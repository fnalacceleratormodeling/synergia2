/*****************************************************************************
*
* Copyright (c) 2000 - 2008, Lawrence Livermore National Security, LLC
* Produced at the Lawrence Livermore National Laboratory
* LLNL-CODE-400142
* All rights reserved.
*
* This file is  part of VisIt. For  details, see https://visit.llnl.gov/.  The
* full copyright notice is contained in the file COPYRIGHT located at the root
* of the VisIt distribution or at http://www.llnl.gov/visit/copyright.html.
*
* Redistribution  and  use  in  source  and  binary  forms,  with  or  without
* modification, are permitted provided that the following conditions are met:
*
*  - Redistributions of  source code must  retain the above  copyright notice,
*    this list of conditions and the disclaimer below.
*  - Redistributions in binary form must reproduce the above copyright notice,
*    this  list of  conditions  and  the  disclaimer (as noted below)  in  the
*    documentation and/or other materials provided with the distribution.
*  - Neither the name of  the LLNS/LLNL nor the names of  its contributors may
*    be used to endorse or promote products derived from this software without
*    specific prior written permission.
*
* THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT  HOLDERS AND CONTRIBUTORS "AS IS"
* AND ANY EXPRESS OR  IMPLIED WARRANTIES, INCLUDING,  BUT NOT  LIMITED TO, THE
* IMPLIED WARRANTIES OF MERCHANTABILITY AND  FITNESS FOR A PARTICULAR  PURPOSE
* ARE  DISCLAIMED. IN  NO EVENT  SHALL LAWRENCE  LIVERMORE NATIONAL  SECURITY,
* LLC, THE  U.S.  DEPARTMENT OF  ENERGY  OR  CONTRIBUTORS BE  LIABLE  FOR  ANY
* DIRECT,  INDIRECT,   INCIDENTAL,   SPECIAL,   EXEMPLARY,  OR   CONSEQUENTIAL
* DAMAGES (INCLUDING, BUT NOT  LIMITED TO, PROCUREMENT OF  SUBSTITUTE GOODS OR
* SERVICES; LOSS OF  USE, DATA, OR PROFITS; OR  BUSINESS INTERRUPTION) HOWEVER
* CAUSED  AND  ON  ANY  THEORY  OF  LIABILITY,  WHETHER  IN  CONTRACT,  STRICT
* LIABILITY, OR TORT  (INCLUDING NEGLIGENCE OR OTHERWISE)  ARISING IN ANY  WAY
* OUT OF THE  USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
* DAMAGE.
*
*****************************************************************************/

// ************************************************************************* //
//                            avtCustomHDF5FileFormat.C                           //
// ************************************************************************* //

#include <avtCustomHDF5FileFormat.h>

#include <list>
#include <string>

#include <vtkFloatArray.h>
#include <vtkRectilinearGrid.h>
#include <vtkStructuredGrid.h>
#include <vtkUnstructuredGrid.h>

#include <avtDatabaseMetaData.h>

#include <DBOptionsAttributes.h>
#include <Expression.h>

#include <InvalidVariableException.h>

#include <DebugStream.h>
//
// Define this symbol BEFORE including hdf5.h to indicate the HDF5 code
// in this file uses version 1.6 of the HDF5 API. This is harmless for
// versions of HDF5 before 1.8 and ensures correct compilation with
// version 1.8 and thereafter. When, and if, the HDF5 code in this file
// is explicitly upgraded to the 1.8 API, this symbol should be removed.
#define H5_USE_16_API
#include <hdf5.h>
#include <visit-hdf5.h>

using     std::string;

namespace CustomHDF5ReaderAux
{
    bool readDatasetDimensions(hid_t file, hsize_t dims[3])
    {
        hid_t coordGroup = H5Gopen(file, "Coordinates");
        if (coordGroup > 0)
        {
            const char *coordDataSetNames[3] = { "r", "theta", "z" };
            bool errorOccured = false;
            for (int i = 0; i < 3; ++i)
            {
                hid_t coordDataSet = H5Dopen(coordGroup, coordDataSetNames[i]);
                if (coordDataSet > 0)
                {
                    hid_t coordDataSpace = H5Dget_space(coordDataSet);
                    if (H5Sis_simple(coordDataSpace))
                    {
                        if (H5Sget_simple_extent_ndims(coordDataSpace) == 1)
                        {
                            H5Sget_simple_extent_dims(coordDataSpace, dims+i, NULL);

                        }
                        else
                        {
                            debug1 << "Expected coordinate arrays to have rank one." << std::endl;
                            errorOccured = true;
                            break;
                        }
                    }
                    else
                    {
                        debug1 << "Coordinate data space is not simple." << std::endl;
                        errorOccured = true;
                        break;
                    }
                    H5Dclose(coordDataSet);
                }
                else
                {
                    debug1 << "Could not open coordinate data set " << coordDataSetNames[i] << std::endl;
                    errorOccured = true;
                    break;
                }
            }
            H5Gclose(coordGroup);
            return !errorOccured;
        }
        else
        {
            debug1 << "Could not open group with coordinates." << std::endl;
            return false;
        }

    }

    vtkFloatArray* readFloatDataSet(hid_t file, const char* dataSetName, int expectedRank, hsize_t* dimsOut)
    {
        hid_t dataSet = H5Dopen(file, dataSetName);
        if (dataSet > 0)
        {
            vtkFloatArray *data = 0;
            bool errorOccured = false;
            hid_t dataSpace = H5Dget_space(dataSet);
            if (H5Sis_simple(dataSpace))
            {
                if (H5Sget_simple_extent_ndims(dataSpace) == expectedRank)
                {
                    H5Sget_simple_extent_dims(dataSpace, dimsOut, NULL);

                    hsize_t numElementsInDataSet = 1;
                    for (int i=0; i<expectedRank; ++i) numElementsInDataSet *= dimsOut[i];

                    data = vtkFloatArray::New();
                    data->SetNumberOfTuples(numElementsInDataSet);
                    float *dataArray = (float *)data->GetVoidPointer(0);
                    H5Dread(dataSet, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, dataArray);
                }
                else
                {
                    debug1 << "Expected coordinate arrays to have rank " << expectedRank <<"." << std::endl;
                }
            }
            else
            {
                debug1 << "Coordinate data space is not simple." << std::endl;
            }

            H5Dclose(dataSet);
            return data;
        }
        else
        {
            debug1 << "Couldn't open data set " << dataSetName << std::endl;
            return 0;
        }
        return 0;
    }

    vtkFloatArray* readFloatDataSetColumn(hid_t file, const char* dataSetName, int expectedRank, int column)
    {
        hid_t dataSet = H5Dopen(file, dataSetName);
        if (dataSet > 0)
        {
            vtkFloatArray *data = 0;
            bool errorOccured = false;
            hid_t dataSpace = H5Dget_space(dataSet);
            if (H5Sis_simple(dataSpace))
            {
                if (H5Sget_simple_extent_ndims(dataSpace) == expectedRank)
                {
                    H5Sget_simple_extent_dims(dataSpace, dimsOut, NULL);

                    hsize_t numElementsInDataSet = 1;
                    for (int i=0; i<expectedRank; ++i) numElementsInDataSet *= dimsOut[i];

                    data = vtkFloatArray::New();
                    data->SetNumberOfTuples(numElementsInDataSet);
                    float *dataArray = (float *)data->GetVoidPointer(0);
                    H5Dread(dataSet, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, dataArray);
                }
                else
                {
                    debug1 << "Expected coordinate arrays to have rank " << expectedRank <<"." << std::endl;
                }
            }
            else
            {
                debug1 << "Coordinate data space is not simple." << std::endl;
            }

            H5Dclose(dataSet);
            return data;
        }
        else
        {
            debug1 << "Couldn't open data set " << dataSetName << std::endl;
            return 0;
        }
        return 0;
    }
}


// ****************************************************************************
//  Method: avtCustomHDF5FileFormat constructor
//
//  Programmer: ghweber -- generated by xml2avt
//  Creation:   Thu Dec 18 14:02:33 PST 2008
//
// ****************************************************************************

avtCustomHDF5FileFormat::avtCustomHDF5FileFormat(const char *filename)
    : avtSTSDFileFormat(filename)
{
}


// ****************************************************************************
//  Method: avtCustomHDF5FileFormat::FreeUpResources
//
//  Purpose:
//      When VisIt is done focusing on a particular timestep, it asks that
//      timestep to free up any resources (memory, file descriptors) that
//      it has associated with it.  This method is the mechanism for doing
//      that.
//
//  Programmer: ghweber -- generated by xml2avt
//  Creation:   Thu Dec 18 14:02:33 PST 2008
//
// ****************************************************************************

void
avtCustomHDF5FileFormat::FreeUpResources(void)
{
}


// ****************************************************************************
//  Method: avtCustomHDF5FileFormat::PopulateDatabaseMetaData
//
//  Purpose:
//      This database meta-data object is like a table of contents for the
//      file.  By populating it, you are telling the rest of VisIt what
//      information it can request from you.
//
//  Programmer: ghweber -- generated by xml2avt
//  Creation:   Thu Dec 18 14:02:33 PST 2008
//
// ****************************************************************************

extern "C"  herr_t
add_var(hid_t loc_id, const char *varname, void *opData)
{
  std::list<std::string> *l = (std::list<std::string>*)opData;
  l->push_back(varname);
}


void
avtCustomHDF5FileFormat::PopulateDatabaseMetaData(avtDatabaseMetaData *md)
{
    //
    // CODE TO ADD A MESH
    //
    string meshname = "mesh";
    //
    // AVT_RECTILINEAR_MESH, AVT_CURVILINEAR_MESH, AVT_UNSTRUCTURED_MESH,
    // AVT_POINT_MESH, AVT_SURFACE_MESH, AVT_UNKNOWN_MESH
    avtMeshType mt;
    mt = AVT_RECTILINEAR_MESH;

    int nblocks = 1;  //<-- this must be 1 for STSD
    int block_origin = 0;
    int spatial_dimension = 3;
    int topological_dimension = 3;
    double *extents = NULL;
    //
    // Here's the call that tells the meta-data object that we have a mesh:
    //
    AddMeshToMetaData(md, meshname, mt, extents, nblocks, block_origin,
            spatial_dimension, topological_dimension);

    std::list<std::string> varList;
    hid_t file = H5Fopen(filename, H5F_ACC_RDONLY, H5P_DEFAULT);
    if (file > 0)
    {
        herr_t idx = H5Giterate(file, "/DataValues", 0, add_var, &varList);
    }
    else
    {
        debug1 << "Couldn't open file " << filename << "." << std::endl;
    }
    for (std::list<std::string>::iterator it = varList.begin(); it!= varList.end(); ++it)
        AddScalarVarToMetaData(md, *it, meshname, AVT_NODECENT);

    //
    // CODE TO ADD A VECTOR VARIABLE
    //
    // string mesh_for_this_var = meshname; // ??? -- could be multiple meshes
    // string varname = ...
    // int vector_dim = 2;
    //
    // AVT_NODECENT, AVT_ZONECENT, AVT_UNKNOWN_CENT
    // avtCentering cent = AVT_NODECENT;
    //
    //
    // Here's the call that tells the meta-data object that we have a var:
    //
    // AddVectorVarToMetaData(md, varname, mesh_for_this_var, cent,vector_dim);
    //

    //
    // CODE TO ADD A TENSOR VARIABLE
    //
    // string mesh_for_this_var = meshname; // ??? -- could be multiple meshes
    // string varname = ...
    // int tensor_dim = 9;
    //
    // AVT_NODECENT, AVT_ZONECENT, AVT_UNKNOWN_CENT
    // avtCentering cent = AVT_NODECENT;
    //
    //
    // Here's the call that tells the meta-data object that we have a var:
    //
    // AddTensorVarToMetaData(md, varname, mesh_for_this_var, cent,tensor_dim);
    //

    //
    // CODE TO ADD A MATERIAL
    //
    // string mesh_for_mat = meshname; // ??? -- could be multiple meshes
    // string matname = ...
    // int nmats = ...;
    // vector<string> mnames;
    // for (int i = 0 ; i < nmats ; i++)
    // {
    //     char str[32];
    //     sprintf(str, "mat%d", i);
    //     -- or -- 
    //     strcpy(str, "Aluminum");
    //     mnames.push_back(str);
    // }
    // 
    // Here's the call that tells the meta-data object that we have a mat:
    //
    // AddMaterialToMetaData(md, matname, mesh_for_mat, nmats, mnames);
    //
    //
    // Here's the way to add expressions:
    //Expression momentum_expr;
    //momentum_expr.SetName("momentum");
    //momentum_expr.SetDefinition("{u, v}");
    //momentum_expr.SetType(Expression::VectorMeshVar);
    //md->AddExpression(&momentum_expr);
    //Expression KineticEnergy_expr;
    //KineticEnergy_expr.SetName("KineticEnergy");
    //KineticEnergy_expr.SetDefinition("0.5*(momentum*momentum)/(rho*rho)");
    //KineticEnergy_expr.SetType(Expression::ScalarMeshVar);
    //md->AddExpression(&KineticEnergy_expr);
}


// ****************************************************************************
//  Method: avtCustomHDF5FileFormat::GetMesh
//
//  Purpose:
//      Gets the mesh associated with this file.  The mesh is returned as a
//      derived type of vtkDataSet (ie vtkRectilinearGrid, vtkStructuredGrid,
//      vtkUnstructuredGrid, etc).
//
//  Arguments:
//      meshname    The name of the mesh of interest.  This can be ignored if
//                  there is only one mesh.
//
//  Programmer: ghweber -- generated by xml2avt
//  Creation:   Thu Dec 18 14:02:33 PST 2008
//
// ****************************************************************************

vtkDataSet *
avtCustomHDF5FileFormat::GetMesh(const char *meshname)
{
    int ndims = 3;
    hsize_t dims[3] = {1,1,1};
    vtkFloatArray *coords[3] = {0,0,0};
    // Read the ndims and number of X,Y,Z nodes from file.

    hid_t file = H5Fopen(filename, H5F_ACC_RDONLY, H5P_DEFAULT);
    if (file > 0)
    {
        coords[0] = CustomHDF5ReaderAux::readFloatDataSet(file, "/Coordinates/r", 1, &dims[0]);
        coords[1] = CustomHDF5ReaderAux::readFloatDataSet(file, "/Coordinates/theta", 1, &dims[1]);
        coords[2] = CustomHDF5ReaderAux::readFloatDataSet(file, "/Coordinates/z", 1, &dims[2]);

        //
        // Create the vtkRectilinearGrid object and set its dimensions
        // and coordinates.
        //
        vtkRectilinearGrid *rgrid = vtkRectilinearGrid::New();
        rgrid->SetDimensions(dims[0], dims[1], dims[2]);
        rgrid->SetXCoordinates(coords[0]);
        coords[0]->Delete();
        rgrid->SetYCoordinates(coords[1]);
        coords[1]->Delete();
        rgrid->SetZCoordinates(coords[2]);
        coords[2]->Delete();
        return rgrid;
    }
    else
    {
        debug5 << "Couldn't open file " << filename << std::endl;
        return 0;
    }
}

// ****************************************************************************
//  Method: avtCustomHDF5FileFormat::GetVar
//
//  Purpose:
//      Gets a scalar variable associated with this file.  Although VTK has
//      support for many different types, the best bet is vtkFloatArray, since
//      that is supported everywhere through VisIt.
//
//  Arguments:
//      varname    The name of the variable requested.
//
//  Programmer: ghweber -- generated by xml2avt
//  Creation:   Thu Dec 18 14:02:33 PST 2008
//
// ****************************************************************************

vtkDataArray *
avtCustomHDF5FileFormat::GetVar(const char *varname)
{
    const char dsPrefix[] = "/DataValues/";
    char dsName[strlen(dsPrefix) + strlen(varname) + 1];
    strcpy(dsName, dsPrefix);
    strcat(dsName, varname);
    hsize_t dim;

    hid_t file = H5Fopen(filename, H5F_ACC_RDONLY, H5P_DEFAULT);
    if (file > 0)
    {
        vtkDataArray *res = CustomHDF5ReaderAux::readFloatDataSet(file, dsName, 1, &dim);
        H5Fclose(file);
        if (!res) EXCEPTION1(InvalidVariableException, varname);
        return res;
    }
    else
    {
        debug5 << "Couldn't open file " << filename << std::endl;
        return 0;
    }
    //
    // If you have a file format where variables don't apply (for example a
    // strictly polygonal format like the STL (Stereo Lithography) format,
    // then uncomment the code below.
    //
    // EXCEPTION1(InvalidVariableException, varname);
    //

    //
    // If you do have a scalar variable, here is some code that may be helpful.
    //
    // int ntuples = XXX; // this is the number of entries in the variable.
    // vtkFloatArray *rv = vtkFloatArray::New();
    // rv->SetNumberOfTuples(ntuples);
    // for (int i = 0 ; i < ntuples ; i++)
    // {
    //      rv->SetTuple1(i, VAL);  // you must determine value for ith entry.
    // }
    //
    // return rv;
    //
}


// ****************************************************************************
//  Method: avtCustomHDF5FileFormat::GetVectorVar
//
//  Purpose:
//      Gets a vector variable associated with this file.  Although VTK has
//      support for many different types, the best bet is vtkFloatArray, since
//      that is supported everywhere through VisIt.
//
//  Arguments:
//      varname    The name of the variable requested.
//
//  Programmer: ghweber -- generated by xml2avt
//  Creation:   Thu Dec 18 14:02:33 PST 2008
//
// ****************************************************************************

vtkDataArray *
avtCustomHDF5FileFormat::GetVectorVar(const char *varname)
{
    return 0;
    //
    // If you have a file format where variables don't apply (for example a
    // strictly polygonal format like the STL (Stereo Lithography) format,
    // then uncomment the code below.
    //
    // EXCEPTION1(InvalidVariableException, varname);
    //

    //
    // If you do have a vector variable, here is some code that may be helpful.
    //
    // int ncomps = YYY;  // This is the rank of the vector - typically 2 or 3.
    // int ntuples = XXX; // this is the number of entries in the variable.
    // vtkFloatArray *rv = vtkFloatArray::New();
    // int ucomps = (ncomps == 2 ? 3 : ncomps);
    // rv->SetNumberOfComponents(ucomps);
    // rv->SetNumberOfTuples(ntuples);
    // float *one_entry = new float[ucomps];
    // for (int i = 0 ; i < ntuples ; i++)
    // {
    //      int j;
    //      for (j = 0 ; j < ncomps ; j++)
    //           one_entry[j] = ...
    //      for (j = ncomps ; j < ucomps ; j++)
    //           one_entry[j] = 0.;
    //      rv->SetTuple(i, one_entry); 
    // }
    //
    // delete [] one_entry;
    // return rv;
    //
}

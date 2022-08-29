#include "space_charge_3d_fd_utils.h"

/* --------------------------------------------------------------------- */
/*!
  Initialize sequential vectors on each MPI rank
  \param   lctx - local context
  \param   gctx - global context
  \return  ierr - PetscErrorCode
  */

PetscErrorCode
init_local_vecs(LocalCtx& lctx, GlobalCtx& gctx)
{
  PetscFunctionBeginUser;

  lctx.seqphi_view = karray1d_dev("seqphi", gctx.nsize);
  lctx.seqrho_view = karray1d_dev("seqphi", gctx.nsize);

  PetscCall(gctx.VecCreate_type_WithArray(PETSC_COMM_SELF,
                                          1,
                                          gctx.nsize,
                                          PETSC_DECIDE,
                                          lctx.seqphi_view.data(),
                                          &lctx.seqphi));
  PetscCall(PetscObjectSetName((PetscObject)(lctx.seqphi), "seqphi_on_lctx"));

  PetscCall(gctx.VecCreate_type_WithArray(PETSC_COMM_SELF,
                                          1,
                                          gctx.nsize,
                                          PETSC_DECIDE,
                                          lctx.seqrho_view.data(),
                                          &lctx.seqrho));
  PetscCall(PetscObjectSetName((PetscObject)(lctx.seqrho), "seqrho_on_lctx"));

  lctx.enx = karray1d_dev("enx", gctx.nsize);
  lctx.eny = karray1d_dev("eny", gctx.nsize);
  lctx.enz = karray1d_dev("enz", gctx.nsize);

  PetscFunctionReturn(0);
}

/* --------------------------------------------------------------------- */
/*!
  Initialize solver subcomms
  \param   sctx - subcomm context
  \param   gctx - global context
  \return  ierr - PetscErrorCode
*/
PetscErrorCode
init_solver_subcomms(SubcommCtx& sctx, GlobalCtx& gctx)
{

  PetscFunctionBeginUser;

  PetscCall(PetscSubcommCreate(PETSC_COMM_WORLD, &(sctx.solverpsubcomm)));
  PetscCall(PetscSubcommSetNumber(sctx.solverpsubcomm, gctx.nsubcomms));
  PetscCall(PetscSubcommSetType(sctx.solverpsubcomm, PETSC_SUBCOMM_CONTIGUOUS));

  PetscCall(PetscSubcommGetChild(sctx.solverpsubcomm, &(sctx.solversubcomm)));
  PetscCall(MPIU_Allreduce(&(gctx.global_rank),
                           &(sctx.solversubcommid),
                           1,
                           MPIU_INT,
                           MPI_MAX,
                           sctx.solversubcomm));
  PetscCallMPI(MPI_Comm_rank(sctx.solversubcomm, &(sctx.solver_rank)));
  PetscCallMPI(MPI_Comm_size(sctx.solversubcomm, &(sctx.solver_size)));

  /* collect the subcommids on each MPI rank & remove duplicates */
  gctx.sids.resize(gctx.global_size);
  PetscCallMPI(MPI_Allgather(&sctx.solversubcommid,
                             1,
                             MPI_INT,
                             gctx.sids.data(),
                             1,
                             MPI_INT,
                             PETSC_COMM_WORLD));
  std::sort(gctx.sids.begin(), gctx.sids.end());
  gctx.sids.erase(std::unique(gctx.sids.begin(), gctx.sids.end()),
                  gctx.sids.end());

  if (gctx.debug) {
    PetscCall(
      PetscPrintf(PETSC_COMM_WORLD, "\nsolver-subcomms have been created!\n"));
    PetscCall(PetscSubcommView(sctx.solverpsubcomm, PETSC_VIEWER_STDOUT_WORLD));
    PetscCall(PetscPrintf(PETSC_COMM_WORLD, "\n"));
    PetscCall(PetscSynchronizedPrintf(
      PETSC_COMM_WORLD,
      "size of gct.sids vector on global-rank %d is : %d\n",
      gctx.global_rank,
      static_cast<int>(gctx.sids.size())));
    PetscCall(PetscSynchronizedFlush(PETSC_COMM_WORLD, PETSC_STDOUT));
    PetscCall(PetscPrintf(PETSC_COMM_WORLD, "\n"));
  }

  PetscFunctionReturn(0);
}

/* --------------------------------------------------------------------- */
/*!
  Initialize vectors on each solver-subcommunicator
  \param   sctx - subcomm context
  \param   gctx - global context
  \return  ierr - PetscErrorCode
  */
PetscErrorCode
init_subcomm_vecs(SubcommCtx& sctx, GlobalCtx& gctx)
{
  PetscInt size;

  PetscFunctionBeginUser;
  PetscCall(VecCreate(sctx.solversubcomm, &sctx.phi_subcomm));
  PetscCall(VecSetType(sctx.phi_subcomm, gctx.vectype));
  PetscCall(VecSetSizes(sctx.phi_subcomm, PETSC_DECIDE, gctx.nsize));
  PetscCall(VecSetFromOptions(sctx.phi_subcomm));
  PetscCall(
    PetscObjectSetName((PetscObject)(sctx.phi_subcomm), "phi_subcomm_on_sctx"));

  PetscCall(VecCreate(sctx.solversubcomm, &sctx.rho_subcomm));
  PetscCall(VecSetType(sctx.rho_subcomm, gctx.vectype));
  PetscCall(VecSetSizes(sctx.rho_subcomm, PETSC_DECIDE, gctx.nsize));
  PetscCall(VecSetFromOptions(sctx.rho_subcomm));
  PetscCall(
    PetscObjectSetName((PetscObject)(sctx.rho_subcomm), "rho_subcomm_on_sctx"));

  if (gctx.debug) {
    PetscCall(VecGetSize(sctx.rho_subcomm, &size));
    PetscCall(
      PetscPrintf(PetscObjectComm((PetscObject)sctx.rho_subcomm),
                  "Hi there from solver-subcomm with id %d, the size of "
                  "subcomm vec here is %d\n",
                  sctx.solversubcommid,
                  size));
    PetscCall(PetscPrintf(PETSC_COMM_WORLD, "\n"));
  }

  PetscFunctionReturn(0);
}

/* --------------------------------------------------------------------- */
/*!
  Initialize DMDA and Matrix to solve Poisson Eq on each solver-subcommunicator
  \param   sctx - subcomm context
  \param   gctx - global context
  \return  ierr - PetscErrorCode
  */
PetscErrorCode
init_subcomm_mat(SubcommCtx& sctx, GlobalCtx& gctx)
{
  PetscFunctionBeginUser;

  PetscCall(DMDACreate3d(gctx.bunch_comm,
                         DM_BOUNDARY_NONE,
                         DM_BOUNDARY_NONE,
                         DM_BOUNDARY_NONE,
                         DMDA_STENCIL_STAR,
                         gctx.nsize_x,
                         gctx.nsize_y,
                         gctx.nsize_z,
                         PETSC_DECIDE,
                         PETSC_DECIDE,
                         PETSC_DECIDE,
                         1,
                         1,
                         NULL,
                         NULL,
                         NULL,
                         &sctx.da));
  PetscCall(DMSetMatType(sctx.da, gctx.mattype));
  PetscCall(DMSetVecType(sctx.da, gctx.vectype));
  PetscCall(DMSetFromOptions(sctx.da));
  PetscCall(DMSetUp(sctx.da));

  /* create discretization matrix */
  PetscCall(DMCreateMatrix(sctx.da, &(sctx.A)));

  PetscFunctionReturn(0);
}

/* --------------------------------------------------------------------- */
/*!
  Compute the RHS matrix to solve the (scaled) Poisson's Eq, setup the
  krylov solver associated with the RHS matrix
  \param   sctx - subcomm context
  \param   gctx - global context
  \return  ierr - PetscErrorCode
  */
PetscErrorCode
compute_mat(SubcommCtx& sctx, GlobalCtx& gctx)
{
  PetscFunctionBeginUser;

  DMDALocalInfo info; /* For storing DMDA info */
  PetscInt i, j, k;
  PetscScalar v[7];
  PetscScalar hx, hy, hz;
  PetscScalar hxhydhz, hxdhyhz, dhxhyhz;
  PetscScalar coordsmin[3], coordsmax[3];
  MatStencil row, col[7];

  PetscFunctionBeginUser;
  PetscCall(DMDAGetLocalInfo(sctx.da, &info));
  PetscCall(DMGetBoundingBox(sctx.da, coordsmin, coordsmax));

  hx = (coordsmax[0] - coordsmin[0]) / (PetscReal)(info.mx);
  hy = (coordsmax[1] - coordsmin[1]) / (PetscReal)(info.my);
  hz = (coordsmax[2] - coordsmin[2]) / (PetscReal)(info.mz);

  hxhydhz = (hx * hy) / hz;
  hxdhyhz = (hx * hz) / hy;
  dhxhyhz = (hy * hz) / hx;

  for (k = info.zs; k < info.zs + info.zm; k++) {
    for (j = info.ys; j < info.ys + info.ym; j++) {
      for (i = info.xs; i < info.xs + info.xm; i++) {
        row.i = i;
        row.j = j;
        row.k = k;
        if (i == 0 || j == 0 || k == 0 || i == info.mx - 1 ||
            j == info.my - 1 || k == info.mz - 1) {
          v[0] = 1.0; // on boundary: trivial equation
          PetscCall(
            MatSetValuesStencil(sctx.A, 1, &row, 1, &row, v, INSERT_VALUES));
        } else {
          v[0] = -hxhydhz;
          col[0].i = i;
          col[0].j = j;
          col[0].k = k - 1;

          v[1] = -hxdhyhz;
          col[1].i = i;
          col[1].j = j - 1;
          col[1].k = k;

          v[2] = -dhxhyhz;
          col[2].i = i - 1;
          col[2].j = j;
          col[2].k = k;

          v[3] = 2 * (hxhydhz + hxdhyhz + dhxhyhz);
          col[3].i = row.i;
          col[3].j = row.j;
          col[3].k = row.k;

          v[4] = -dhxhyhz;
          col[4].i = i + 1;
          col[4].j = j;
          col[4].k = k;

          v[5] = -hxdhyhz;
          col[5].i = i;
          col[5].j = j + 1;
          col[5].k = k;

          v[6] = -hxhydhz;
          col[6].i = i;
          col[6].j = j;
          col[6].k = k + 1;

          PetscCall(
            MatSetValuesStencil(sctx.A, 1, &row, 7, col, v, INSERT_VALUES));
        }
      }
    }
  }

  PetscCall(MatAssemblyBegin(sctx.A, MAT_FINAL_ASSEMBLY));
  PetscCall(MatAssemblyEnd(sctx.A, MAT_FINAL_ASSEMBLY));

  /* create krylov solver */
  PetscCall(KSPCreate(sctx.solversubcomm, &sctx.ksp));
  PetscCall(KSPSetType(sctx.ksp, KSPGMRES));
  PetscCall(KSPSetOperators(sctx.ksp, sctx.A, sctx.A));
  PetscCall(KSPSetFromOptions(sctx.ksp));
  PetscCall(KSPSetUp(sctx.ksp));
  // PetscCall(KSPSetDM(sctx.ksp, sctx.da));

  /* set preconditioner */
  PetscCall(KSPGetPC(sctx.ksp, &(sctx.pc)));
  PetscCall(PCSetType(sctx.pc, PCASM));
  PetscCall(PCSetUp(sctx.pc));

  /* Enable KSP logging if options are set */
  if (gctx.ksp_view) {
    PetscCall(KSPView(
      sctx.ksp, PETSC_VIEWER_STDOUT_(PetscObjectComm((PetscObject)sctx.ksp))));
  }
  if (gctx.ksp_converged_reason) {
    PetscCall(KSPConvergedReasonView(
      sctx.ksp, PETSC_VIEWER_STDOUT_(PetscObjectComm((PetscObject)sctx.ksp))));
  }

  if (gctx.ksp_monitor_residual) {
    PetscCall(KSPMonitorSet(sctx.ksp, &(MyMonitor), PETSC_NULL, PETSC_NULL));
  }

  PetscFunctionReturn(0);
}
/* --------------------------------------------------------------------- */
/*!
  Solve the scaled Poisson Eq!
  \param   sctx - subcomm context
  \param   gctx - global context
  \return  ierr - PetscErrorCode
  */
PetscErrorCode
solve(SubcommCtx& sctx, GlobalCtx& gctx)
{

  DMDALocalInfo info; /* For storing DMDA info */
  PetscScalar hx, hy, hz;
  PetscScalar coordsmin[3], coordsmax[3];

  PetscFunctionBeginUser;
  PetscCall(DMDAGetLocalInfo(sctx.da, &info));
  PetscCall(DMGetBoundingBox(sctx.da, coordsmin, coordsmax));

  hx = (coordsmax[0] - coordsmin[0]) / (PetscReal)(info.mx);
  hy = (coordsmax[1] - coordsmin[1]) / (PetscReal)(info.my);
  hz = (coordsmax[2] - coordsmin[2]) / (PetscReal)(info.mz);

  /* Scaling factor of hx*hy*hz */
  PetscCall(VecScale(sctx.rho_subcomm, hx * hy * hz));

  /* Solve for phi! */
  PetscCall(KSPSolve(sctx.ksp, sctx.rho_subcomm, sctx.phi_subcomm));

  PetscFunctionReturn(0);
}

/* --------------------------------------------------------------------- */
/*!
  Initialize global (alias of local) to subcomm scatters
  \param   sctx - subcomm context
  \param   gctx - global context
  \return  ierr - PetscErrorCode
*/
PetscErrorCode
init_global_subcomm_scatters(SubcommCtx& sctx, GlobalCtx& gctx)
{

  PetscInt start;
  PetscInt localsize;

  PetscFunctionBeginUser;

  /* localsize should be gctx->nsize, but adding a check doesn't hurt */
  PetscCall(VecGetLocalSize(gctx.rho_global_local, &localsize));
  PetscCall(VecGetOwnershipRange(gctx.rho_global_local, &start, NULL));

  PetscCall(ISCreateStride(
    PETSC_COMM_SELF, localsize, start, 1, &gctx.ix_scat_glocal_to_subcomms));
  PetscCall(ISCreateStride(
    PETSC_COMM_SELF, localsize, 0, 1, &gctx.iy_scat_glocal_to_subcomms));

  PetscCall(PetscObjectSetName((PetscObject)(gctx.ix_scat_glocal_to_subcomms),
                               "ix_scat_glocal_to_subcomms"));
  PetscCall(PetscObjectSetName((PetscObject)(gctx.iy_scat_glocal_to_subcomms),
                               "iy_scat_glocal_to_subcomms"));

  /* resize to ensure we only create as many scatters
     as the number of subcomms ! */
  gctx.scat_glocal_to_subcomms.resize(gctx.nsubcomms);

  for (PetscInt i = 0; i < gctx.nsubcomms; i++) {
    PetscCall(VecScatterCreate(gctx.rho_global_local,
                               gctx.ix_scat_glocal_to_subcomms,
                               gctx.rho_global_subcomm[i],
                               gctx.iy_scat_glocal_to_subcomms,
                               &(gctx.scat_glocal_to_subcomms[i])));
  }

  if (gctx.debug) {
    for (PetscInt i = 0; i < gctx.nsubcomms; i++) {
      PetscCall(PetscPrintf(PETSC_COMM_WORLD, "global-to-subcomm scatter\n"));
      PetscCall(VecScatterView(gctx.scat_glocal_to_subcomms[i],
                               PETSC_VIEWER_STDOUT_WORLD));
    }
    PetscCall(PetscPrintf(PETSC_COMM_WORLD, "\n"));
  }

  PetscFunctionReturn(0);
}

/* --------------------------------------------------------------------- */
/*!
  Initialize subcomm (alias of local) to local scatters
  \param   lctx - local context
  \param   sctx - subcomm context
  \param   gctx - global context
  \return  ierr - PetscErrorCode
*/
PetscErrorCode
init_subcomm_local_scatters(LocalCtx& lctx, SubcommCtx& sctx, GlobalCtx& gctx)
{

  PetscInt start;
  PetscInt localsize;

  PetscFunctionBeginUser;

  /* localsize should be gctx->nsize, but adding a check doesn't hurt */
  PetscCall(VecGetLocalSize(sctx.phi_subcomm_local, &localsize));
  PetscCall(VecGetOwnershipRange(sctx.phi_subcomm_local, &start, NULL));

  PetscCall(ISCreateStride(
    PETSC_COMM_SELF, localsize, start, 1, &sctx.ix_scat_subcomms_to_local));
  PetscCall(ISCreateStride(
    PETSC_COMM_SELF, localsize, 0, 1, &sctx.iy_scat_subcomms_to_local));

  /* This scatter will always be run as a SCATTER_REVERSE,
     perhaps the naming terminilogy may be updated to prevent
     any confusion in the future.*/
  PetscCall(VecScatterCreate(sctx.phi_subcomm_local,
                             sctx.ix_scat_subcomms_to_local,
                             sctx.phi_subcomm,
                             sctx.iy_scat_subcomms_to_local,
                             &sctx.scat_subcomm_to_local));

  if (gctx.debug) {
    PetscCall(PetscPrintf(PETSC_COMM_WORLD, "subcomm-to-local scatter\n"));
    PetscCall(VecScatterView(
      sctx.scat_subcomm_to_local,
      PETSC_VIEWER_STDOUT_(PetscObjectComm((PetscObject)sctx.phi_subcomm))));
    PetscCall(PetscPrintf(PETSC_COMM_WORLD, "\n"));
    PetscCall(PetscBarrier(NULL));
  }

  PetscFunctionReturn(0);
}

/* --------------------------------------------------------------------- */
/*!
  finalize by destroying data structures
  \param   lctx - local context
  \param   sctx - subcomm context
  \param   gctx - global context
  \return  ierr - PetscErrorCode
*/
PetscErrorCode
finalize(LocalCtx& lctx, SubcommCtx& sctx, GlobalCtx& gctx)
{

  PetscFunctionBeginUser;

  /* Destroy global aliases of local vectors */
  PetscCall(VecDestroy(&(gctx.phi_global_local)));
  PetscCall(VecDestroy(&(gctx.rho_global_local)));

  /* Destroy global aliases of subcomm vectors */
  for (PetscInt i = 0; i < gctx.nsubcomms; i++) {
    PetscCall(VecDestroy(&(gctx.phi_global_subcomm[i])));
    PetscCall(VecDestroy(&(gctx.rho_global_subcomm[i])));
  }

  /* Destroy subcomm vectors */
  PetscCall(VecDestroy(&(sctx.phi_subcomm)));
  PetscCall(VecDestroy(&(sctx.rho_subcomm)));

  /* Destroy global (alias of local) to subcomm scatters */
  for (PetscInt i = 0; i < gctx.nsubcomms; i++) {
    PetscCall(VecScatterDestroy(&(gctx.scat_glocal_to_subcomms[i])));
  }
  PetscCall(ISDestroy(&gctx.ix_scat_glocal_to_subcomms));
  PetscCall(ISDestroy(&gctx.iy_scat_glocal_to_subcomms));

  /* Destroy subcomm aliases of local vectors */
  PetscCall(VecDestroy(&sctx.phi_subcomm_local));
  PetscCall(VecDestroy(&sctx.rho_subcomm_local));

  /* Destroy DMDA and matrix on subcomm */
  PetscCall(KSPDestroy(&sctx.ksp));
  PetscCall(MatDestroy(&sctx.A));
  PetscCall(DMDestroy(&sctx.da));

  /* Destroy subcomm vectors */
  PetscCall(VecDestroy(&(lctx.seqphi)));
  PetscCall(VecDestroy(&(lctx.seqrho)));

  /* Destroy subcomm (alias of local) to local scatters */
  PetscCall(VecScatterDestroy(&sctx.scat_subcomm_to_local));
  PetscCall(ISDestroy(&sctx.ix_scat_subcomms_to_local));
  PetscCall(ISDestroy(&sctx.iy_scat_subcomms_to_local));

  /* Destroy subcomms */
  PetscCall(PetscSubcommDestroy(&(sctx.solverpsubcomm)));
  PetscFunctionReturn(0);
}

/* --------------------------------------------------------------------- */
/*!
  MyMonitor wrapper around KSPMonitorResidual,
  see https://petsc.org/main/docs/manualpages/KSP/KSPMonitorSet/
  \param   KSP - ksp solver object
  \param   n - iteration number
  \param   rnorm - (estimated) 2-norm of (preconditioned) residual
  \param   mctx - optional monitoring context, as set by KSPMonitorSet()
  \return  ierr - PetscErrorCode
*/
PetscErrorCode
MyMonitor(KSP ksp, PetscInt it, PetscReal rnorm, void* mctx)
{
  PetscFunctionBeginUser;
  PetscViewerAndFormat* vf;
  PetscCall(PetscViewerAndFormatCreate(
    PETSC_VIEWER_STDOUT_(PetscObjectComm((PetscObject)ksp)),
    PETSC_VIEWER_ASCII_INFO,
    &vf));
  PetscCall(KSPMonitorResidual(ksp, it, rnorm, vf));
  PetscFunctionReturn(0);
}

template<class T>
Array_3d<T>::Array_3d() : Array_nd<T>()
{
    this->shape = vector3(0,0,0);
    this->strides = vector3(0,0,0);
}

template<class T>
Array_3d<T>::Array_3d(const int nx, const int ny, const int nz, T *data_ptr) :
    Array_nd<T>()
{
    this->construct(vector3(nx,ny,nz),data_ptr);
}

template<class T>
Array_3d<T>::Array_3d(const int nx, const int ny, const int nz,
        const int stride_x, const int stride_y, const int stride_z,
        T *data_ptr) :
    Array_nd<T>()
{
    this->construct(vector3(nx,ny,nz),vector3(stride_x,stride_y,stride_z),
        data_ptr);
}

template<class T>
Array_3d<T>::Array_3d(const Array_nd<T>& original)
{
    //~ std::cout << "converting from nd to 3d\n";
    if (original.get_rank() != 3) {       
        throw
            std::runtime_error("Attempt to convert Array_nd of rank !=3 to Array_3d");
    }
    copy_construct(original);
}

template<class T>
void
Array_3d<T>::reshape(const int nx, const int ny, const int nz, T *data_ptr)
{
    bool shape_changed = this->different_shape(vector3(nx,ny,nz));
    if (shape_changed && this->shape_frozen) {
        throw
        std::out_of_range("Attempt to change the shape of a frozen Array_3d");
    }
    if (shape_changed) {
        this->construct(vector3(nx,ny,nz),data_ptr);
    } else {
    	if (data_ptr != 0) {
    		this->data_ptr = data_ptr;
    	}
    }
}

template<class T>
void
Array_3d<T>::reshape(const int nx, const int ny, const int nz,
        const int stride_x, const int stride_y, const int stride_z, T *data_ptr)
{
    bool shape_changed = this->different_shape(vector3(nx,ny,nz));
    if (shape_changed && this->shape_frozen) {
        throw
        std::out_of_range("Attempt to change the shape of a frozen Array_3d");
    }
    if (shape_changed) {
        this->construct(vector3(nx,ny,nz),vector3(stride_x,stride_y,stride_z),data_ptr);
    } else {
    	if (data_ptr != 0) {
    		this->data_ptr = data_ptr;
    	}
    }
}

template<class T>
inline int
Array_3d<T>::offset(const int i, const int j, const int k) const
{
    return i*this->strides[0] + j*this->strides[1] + k*this->strides[2];
}


template<class T>
inline T&
Array_3d<T>::at(const int i, const int j, const int k)
{
    if (bounds_check(i,j,k)) {
        return this->data_ptr[offset(i,j,k)];
    } else {
        throw
        std::out_of_range("Bounds check failed");
    }
}

template<class T>
inline T
Array_3d<T>::at(const int i, const int j, const int k) const
{
    if (bounds_check(i,j,k)) {
        return this->data_ptr[offset(i,j,k)];
    } else {
        throw
        std::out_of_range("Bounds check failed");
    }
}

template<class T>
inline T&
Array_3d<T>::operator()(const int i, const int j, const int k)
{
    return this->data_ptr[offset(i,j,k)];
}

template<class T>
inline T
Array_3d<T>::operator()(const int i, const int j, const int k) const
{
    return this->data_ptr[offset(i,j,k)];
}

template<class T>
inline bool
Array_3d<T>::bounds_check(const int i, const int j, const int k) const
{
    if ( (i<0) || (i>=this->shape[0]) ||
        (j<0) || (j>=this->shape[1]) ||
        (k<0) || (k>=this->shape[2]) ) {
        return false;
    }
    return true;
}

template<class T>
Array_nd<T>
Array_3d<T>::slice(Range i_range, Range j_range, Range k_range)
{
	return Array_nd<T>::slice(vector3(i_range,j_range,k_range));
}

template<class T>
Array_3d<T>::~Array_3d()
{
}

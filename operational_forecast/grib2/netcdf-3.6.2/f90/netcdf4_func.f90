  !
  ! NetCDF-4 extra routines:
  !
  ! -----------
  function nf90_create_par(path, cmode, comm, info, ncid)
    character (len = *), intent(in) :: path
    integer, intent(in) :: cmode
    integer, intent(in) :: comm
    integer, intent(in) :: info
    integer, intent(out) :: ncid
    integer :: nf90_create_par
    nf90_create_par = nf_create_par(path, cmode, comm, info, ncid)
  end function nf90_create_par
  ! -----------
  function nf90_open_par(path, cmode, comm, info, ncid)
    character (len = *), intent(in) :: path
    integer, intent(in) :: cmode
    integer, intent(in) :: comm
    integer, intent(in) :: info
    integer, intent(out) :: ncid
    integer :: nf90_open_par
  
    nf90_open_par = nf_open_par(path, cmode, comm, info, ncid)
  end function nf90_open_par
  ! -----------
  function nf90_var_par_access(ncid, varid, access)
    integer, intent(in) :: ncid
    integer, intent(in) :: varid
    integer, intent(in) :: access
    integer :: nf90_var_par_access
  
    nf90_var_par_access = nf_var_par_access(ncid, varid, access)
  end function nf90_var_par_access
  ! -----------
  function nf90_inq_ncid(ncid, name, grp_ncid)
    integer, intent(in) :: ncid
    character (len = *), intent(in) :: name
    integer, intent(out) :: grp_ncid
    integer :: nf90_inq_ncid
  
    nf90_inq_ncid = nf_inq_ncid(ncid, name, grp_ncid)
  end function nf90_inq_ncid
  ! -----------
  function nf90_inq_grps(ncid, numgrps, ncids)
    integer, intent(in) :: ncid
    integer, intent(out) :: numgrps
    integer, intent(out) :: ncids
    integer :: nf90_inq_grps
  
    nf90_inq_grps = nf_inq_grps(ncid, numgrps, ncids)
  end function nf90_inq_grps
  ! -----------
  function nf90_inq_grpname(ncid, name)
    integer, intent(in) :: ncid
    character (len = *), intent(out) :: name
    integer :: nf90_inq_grpname
  
    nf90_inq_grpname = nf_inq_grpname(ncid, name)
  end function nf90_inq_grpname
  ! -----------
  function nf90_inq_varids(ncid, nvars, varids)
    integer, intent(in) :: ncid
    integer, intent(out) :: nvars
    integer, intent(out) :: varids
    integer :: nf90_inq_varids
  
    nf90_inq_varids = nf_inq_varids(ncid, nvars, varids)
  end function nf90_inq_varids
  ! -----------
  function nf90_inq_dimids(ncid, ndims, dimids, include_parents)
    integer, intent(in) :: ncid
    integer, intent(out) :: ndims
    integer, intent(out) :: dimids
    integer, intent(out) :: include_parents
    integer :: nf90_inq_dimids
  
    nf90_inq_dimids = nf_inq_dimids(ncid, ndims, dimids, include_parents)
  end function nf90_inq_dimids
  ! -----------
  function nf90_inq_typeids(ncid, ntypes, typeids)
    integer, intent(in) :: ncid
    integer, intent(out) :: ntypes
    integer, intent(out) :: typeids
    integer :: nf90_inq_typeids
  
    nf90_inq_typeids = nf_inq_typeids(ncid, ntypes, typeids)
  end function nf90_inq_typeids
  ! -----------
  function nf90_def_grp(parent_ncid, name, new_ncid)
    integer, intent(in) :: parent_ncid
    character (len = *), intent(in) :: name
    integer, intent(out) :: new_ncid
    integer :: nf90_def_grp
  
    nf90_def_grp = nf_def_grp(parent_ncid, name, new_ncid)
  end function nf90_def_grp
  ! -----------
  function nf90_def_compound(ncid, size, name, typeid)
    integer, intent(in) :: ncid
    integer, intent(in) :: size
    character (len = *), intent(in) :: name
    integer, intent(out) :: typeid
    integer :: nf90_def_compound
  
    nf90_def_compound = nf_def_compound(ncid, size, name, typeid)
  end function nf90_def_compound
  ! -----------
  function nf90_insert_compound(ncid, xtype, name, offset, field_typeid)
    integer, intent(in) :: ncid
    integer, intent(in) :: xtype
    character (len = *), intent(in) :: name
    integer, intent(in) :: offset
    integer, intent(in) :: field_typeid
    integer :: nf90_insert_compound
  
    nf90_insert_compound = nf_insert_compound(ncid, xtype, name, offset, field_typeid)
  end function nf90_insert_compound
  ! -----------
  function nf90_insert_array_compound(ncid, xtype, name, offset, field_typeid, &
       ndims, dim_sizes)
    integer, intent(in) :: ncid
    integer, intent(in) :: xtype
    character (len = *), intent(in) :: name
    integer, intent(in) :: offset
    integer, intent(in) :: field_typeid
    integer, intent(in) :: ndims
    integer, intent(in) :: dim_sizes
    integer :: nf90_insert_array_compound
  
    nf90_insert_array_compound = nf_insert_array_compound(ncid, xtype, name, &
         offset, field_typeid, ndims, dim_sizes)
  end function nf90_insert_array_compound
  ! -----------
  function nf90_inq_type(ncid, xtype, name, size, nfields)
    integer, intent(in) :: ncid
    integer, intent(in) :: xtype
    character (len = *), intent(out) :: name
    integer, intent(out) :: size
    integer, intent(out) :: nfields
    integer :: nf90_inq_type
  
    nf90_inq_type = nf_inq_type(ncid, xtype, name, size, nfields)
  end function nf90_inq_type
  ! -----------
  function nf90_inq_compound(ncid, xtype, name, size, nfields)
    integer, intent(in) :: ncid
    integer, intent(in) :: xtype
    character (len = *), intent(out) :: name
    integer, intent(out) :: size
    integer, intent(out) :: nfields
    integer :: nf90_inq_compound
  
    nf90_inq_compound = nf_inq_compound(ncid, xtype, name, size, nfields)
  end function nf90_inq_compound
!   ! -----------
  function nf90_inq_compound_field(ncid, xtype, fieldid, name, offset, &
       field_typeid, ndims, dim_sizes)
    integer, intent(in) :: ncid
    integer, intent(in) :: xtype
    integer, intent(in) :: fieldid
    character (len = *), intent(out) :: name
    integer, intent(out) :: offset
    integer, intent(out) :: field_typeid
    integer, intent(out) :: ndims
    integer, intent(out) :: dim_sizes
    integer :: nf90_inq_compound_field
  
    nf90_inq_compound_field = nf_inq_compound_field(ncid, xtype, fieldid, name, offset, &
       field_typeid, ndims, dim_sizes)
  end function nf90_inq_compound_field
  ! -----------
  function nf90_def_vlen(ncid, name, base_typeid, xtypeid)
    integer, intent(in) :: ncid
    character (len = *), intent(in) :: name
    integer, intent(in) :: base_typeid
    integer, intent(out) :: xtypeid
    integer :: nf90_def_vlen
  
    nf90_def_vlen = nf_def_vlen(ncid, name, base_typeid, xtypeid)
  end function nf90_def_vlen
  ! -----------
  function nf90_inq_vlen(ncid, xtype, name, datum_size, base_nc_type)
    integer, intent(in) :: ncid
    integer, intent(in) :: xtype
    character (len = *), intent(out) :: name
    integer, intent(out) :: datum_size
    integer, intent(out) :: base_nc_type
    integer :: nf90_inq_vlen
  
    nf90_inq_vlen = nf_inq_vlen(ncid, xtype, name, datum_size, base_nc_type)
  end function nf90_inq_vlen
!   ! -----------
  function nf90_def_enum(ncid, base_typeid, name, typeid)
    integer, intent(in) :: ncid
    integer, intent(in) :: base_typeid
    character (len = *), intent(in) :: name
    integer, intent(out) :: typeid
    integer :: nf90_def_enum
  
    nf90_def_enum = nf_def_enum(ncid, base_typeid, name, typeid)
  end function nf90_def_enum
  ! -----------
  function nf90_insert_enum(ncid, xtype, name, value)
    integer, intent(in) :: ncid
    integer, intent(in) :: xtype
    character (len = *), intent(in) :: name
    integer, intent(in) :: value
    integer :: nf90_insert_enum
  
    nf90_insert_enum = nf_insert_enum(ncid, xtype, name, value)
  end function nf90_insert_enum
  ! -----------
  function nf90_inq_enum(ncid, xtype, name, base_nc_type, base_size, num_members)
    integer, intent(in) :: ncid
    integer, intent(in) :: xtype
    character (len = *), intent(out) :: name
    integer, intent(out) :: base_nc_type
    integer, intent(out) :: base_size
    integer, intent(out) :: num_members
    integer :: nf90_inq_enum
  
    nf90_inq_enum = nf_inq_enum(ncid, xtype, name, base_nc_type, base_size, num_members)
  end function nf90_inq_enum
  ! -----------
  function nf90_inq_enum_member(ncid, xtype, idx, name, value)
    integer, intent(in) :: ncid
    integer, intent(in) :: xtype
    integer, intent(in) :: idx
    character (len = *), intent(out) :: name
    integer, intent(in) :: value
    integer :: nf90_inq_enum_member
  
    nf90_inq_enum_member = nf_inq_enum_member(ncid, xtype, idx, name, value)
  end function nf90_inq_enum_member
  ! -----------
  function nf90_def_opaque(ncid, size, name, xtype)
    integer, intent(in) :: ncid
    integer, intent(in) :: size
    character (len = *), intent(in) :: name
    integer, intent(out) :: xtype
    integer :: nf90_def_opaque
  
    nf90_def_opaque = nf_def_opaque(ncid, size, name, xtype)
  end function nf90_def_opaque
  ! -----------
  function nf90_inq_opaque(ncid, xtype, name, size)
    integer, intent(in) :: ncid
    integer, intent(in) :: xtype
    character (len = *), intent(out) :: name
    integer, intent(out) :: size
    integer :: nf90_inq_opaque
  
    nf90_inq_opaque = nf_inq_opaque(ncid, xtype, name, size)
  end function nf90_inq_opaque
  ! -----------
  function nf90_def_var_chunking(ncid, varid, chunkalg, chunksizes, extend_increments)
    integer, intent(in) :: ncid
    integer, intent(in) :: varid
    integer, intent(in) :: chunkalg
    integer, intent(in) :: chunksizes
    integer, intent(in) :: extend_increments
    integer :: nf90_def_var_chunking
  
    nf90_def_var_chunking = nf_def_var_chunking(ncid, varid, chunkalg, chunksizes, extend_increments)
  end function nf90_def_var_chunking
  ! -----------
  function nf90_def_var_deflate(ncid, varid, shuffle, deflate, deflate_level)
    integer, intent(in) :: ncid
    integer, intent(in) :: varid
    integer, intent(in) :: shuffle
    integer, intent(in) :: deflate
    integer, intent(in) :: deflate_level
    integer :: nf90_def_var_deflate
  
    nf90_def_var_deflate = nf_def_var_deflate(ncid, varid, shuffle, deflate, deflate_level)
  end function nf90_def_var_deflate
  ! -----------
  function nf90_def_var_fletcher32(ncid, varid, fletcher32)
    integer, intent(in) :: ncid
    integer, intent(in) :: varid
    integer, intent(in) :: fletcher32
    integer :: nf90_def_var_fletcher32
  
    nf90_def_var_fletcher32 = nf_def_var_fletcher32(ncid, varid, fletcher32)
  end function nf90_def_var_fletcher32
  ! -----------
  function nf90_inq_var_chunking(ncid, varid, chunkalg, chunksizes, extend_increments)
    integer, intent(in) :: ncid
    integer, intent(in) :: varid
    integer, intent(out) :: chunkalg
    integer, intent(out) :: chunksizes
    integer, intent(out) :: extend_increments
    integer :: nf90_inq_var_chunking
  
    nf90_inq_var_chunking = nf_inq_var_chunking(ncid, varid, chunkalg, chunksizes, extend_increments)
  end function nf90_inq_var_chunking
  ! -----------
  function nf90_inq_var_deflate(ncid, varid, shuffle, deflate, deflate_level)
    integer, intent(in) :: ncid
    integer, intent(in) :: varid
    integer, intent(out) :: shuffle
    integer, intent(out) :: deflate
    integer, intent(out) :: deflate_level
    integer :: nf90_inq_var_deflate
  
    nf90_inq_var_deflate = nf_inq_var_deflate(ncid, varid, shuffle, deflate, deflate_level)
  end function nf90_inq_var_deflate
  ! -----------
  function nf90_inq_var_fletcher32(ncid, varid, fletcher32)
    integer, intent(in) :: ncid
    integer, intent(in) :: varid
    integer, intent(out) :: fletcher32
    integer :: nf90_inq_var_fletcher32
  
    nf90_inq_var_fletcher32 = nf_inq_var_fletcher32(ncid, varid, fletcher32)
  end function nf90_inq_var_fletcher32
  ! -----------
  function nf90_def_var_endian(ncid, varid, endian)
    integer, intent(in) :: ncid
    integer, intent(in) :: varid
    integer, intent(in) :: endian
    integer :: nf90_def_var_endian
  
    nf90_def_var_endian = nf_def_var_endian(ncid, varid, endian)
  end function nf90_def_var_endian
  ! -----------
  function nf90_inq_var_endian(ncid, varid, endian)
    integer, intent(in) :: ncid
    integer, intent(in) :: varid
    integer, intent(out) :: endian
    integer :: nf90_inq_var_endian
  
    nf90_inq_var_endian = nf_inq_var_endian(ncid, varid, endian)
  end function nf90_inq_var_endian
  ! -----------

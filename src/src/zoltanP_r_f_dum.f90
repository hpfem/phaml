! dummy version of routines in zoltanParams_read_file.  To be used
! with Zoltan Versions <3.0

subroutine Zf90_zoltanparams_read_file(lb_addr,nbytes,ifile,nchar, &
                                       communicator,ierr)
integer, intent(in) :: lb_addr(*), ifile(*)
integer, intent(in) :: nbytes, nchar, communicator
integer, intent(out) :: ierr
ierr=0
end subroutine Zf90_zoltanparams_read_file

subroutine zoltanparams_using_graph(retval)
integer, intent(out) :: retval
retval=0
end subroutine zoltanparams_using_graph

subroutine Zoltan_Get_Struct_Addr(lb,lb_addr)
integer :: lb(*),lb_addr(*)
end subroutine Zoltan_Get_Struct_Addr

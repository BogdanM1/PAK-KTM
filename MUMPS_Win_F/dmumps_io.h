/*

   THIS FILE IS PART OF MUMPS VERSION 4.6.2
   This Version was built on Fri Apr 14 14:59:20 2006


  This version of MUMPS is provided to you free of charge. It is public
  domain, based on public domain software developed during the Esprit IV
  European project PARASOL (1996-1999) by CERFACS, ENSEEIHT-IRIT and RAL. 
  Since this first public domain version in 1999, the developments are
  supported by the following institutions: CERFACS, ENSEEIHT-IRIT, and
  INRIA.

  Main contributors are Patrick Amestoy, Iain Duff, Abdou Guermouche,
  Jacko Koster, Jean-Yves L'Excellent, and Stephane Pralet.

  Up-to-date copies of the MUMPS package can be obtained
  from the Web pages http://www.enseeiht.fr/apo/MUMPS/
  or http://graal.ens-lyon.fr/MUMPS


   THIS MATERIAL IS PROVIDED AS IS, WITH ABSOLUTELY NO WARRANTY
   EXPRESSED OR IMPLIED.  ANY USE IS AT YOUR OWN RISK.


  User documentation of any code that uses this software can
  include this complete notice. You can acknowledge (using
  references [1], [2], and [3] the contribution of this package
  in any scientific publication dependent upon the use of the
  package. You shall use reasonable endeavours to notify
  the authors of the package of this publication.

   [1] P. R. Amestoy, I. S. Duff and  J.-Y. L'Excellent (1998),
   Multifrontal parallel distributed symmetric and unsymmetric solvers,
   in Comput. Methods in Appl. Mech. Eng., 184,  501-520 (2000).

   [2] P. R. Amestoy, I. S. Duff, J. Koster and  J.-Y. L'Excellent,
   A fully asynchronous multifrontal solver using distributed dynamic
   scheduling, SIAM Journal of Matrix Analysis and Applications,
   Vol 23, No 1, pp 15-41 (2001).

   [3] P. R. Amestoy and A. Guermouche and J.-Y. L'Excellent and
   S. Pralet, Hybrid scheduling for the parallel solution of linear
   systems. Parallel Computing Vol 32 (2), pp 136-156 (2006).

*/
/*    $Id: dmumps_io.h,v 1.20 2006/03/27 16:46:59 jylexcel Exp $ */

#if defined(_WIN32) || defined (UPPER)
#define dmumps_is_there_finished_request DMUMPS_IS_THERE_FINISHED_REQUEST
#define dmumps_clean_request DMUMPS_CLEAN_REQUEST
#define dmumps_test_request DMUMPS_TEST_REQUEST
#define dmumps_wait_request DMUMPS_WAIT_REQUEST
#define dmumps_wait_all_requests DMUMPS_WAIT_ALL_REQUESTS
#define dmumps_low_level_init_ooc_c DMUMPS_LOW_LEVEL_INIT_OOC_C
#define dmumps_low_level_write_ooc_c DMUMPS_LOW_LEVEL_WRITE_OOC_C
#define dmumps_low_level_read_ooc_c DMUMPS_LOW_LEVEL_READ_OOC_C
#define dmumps_low_level_direct_read DMUMPS_LOW_LEVEL_DIRECT_READ
#define dmumps_clean_io_data_c DMUMPS_CLEAN_IO_DATA_C
#define dmumps_get_max_nb_req DMUMPS_GET_MAX_NB_REQ
#define dmumps_get_max_file_size DMUMPS_GET_MAX_FILE_SIZE
#define dmumps_ooc_get_nb_files DMUMPS_OOC_GET_NB_FILES
#define dmumps_ooc_get_file_name DMUMPS_OOC_GET_FILE_NAME
#define dmumps_ooc_set_file_name DMUMPS_OOC_SET_FILE_NAME
#define dmumps_ooc_start_low_level DMUMPS_OOC_START_LOW_LEVEL
#define dmumps_ooc_alloc_pointers DMUMPS_OOC_ALLOC_POINTERS
#define dmumps_ooc_print_stats DMUMPS_OOC_PRINT_STATS
#define dmumps_ooc_remove_file DMUMPS_OOC_REMOVE_FILE
#define dmumps_ooc_init_vars DMUMPS_OOC_INIT_VARS
#define dmumps_ooc_end_write DMUMPS_OOC_END_WRITE
#elif defined(Add__)
#define dmumps_is_there_finished_request dmumps_is_there_finished_request__
#define dmumps_clean_request dmumps_clean_request__
#define dmumps_test_request dmumps_test_request__
#define dmumps_wait_request dmumps_wait_request__
#define dmumps_wait_all_requests dmumps_wait_all_requests__
#define dmumps_low_level_init_ooc_c dmumps_low_level_init_ooc_c__
#define dmumps_low_level_write_ooc_c dmumps_low_level_write_ooc_c__
#define dmumps_low_level_read_ooc_c dmumps_low_level_read_ooc_c__
#define dmumps_low_level_direct_read dmumps_low_level_direct_read__
#define dmumps_clean_io_data_c dmumps_clean_io_data_c__
#define dmumps_get_max_nb_req dmumps_get_max_nb_req__
#define dmumps_get_max_file_size dmumps_get_max_file_size__ 
#define dmumps_ooc_get_nb_files dmumps_ooc_get_nb_files__
#define dmumps_ooc_get_file_name dmumps_ooc_get_file_name__
#define dmumps_ooc_set_file_name dmumps_ooc_set_file_name__
#define dmumps_ooc_start_low_level dmumps_ooc_start_low_level__
#define dmumps_ooc_alloc_pointers dmumps_ooc_alloc_pointers__
#define dmumps_ooc_print_stats dmumps_ooc_print_stats__
#define dmumps_ooc_remove_file dmumps_ooc_remove_file__
#define dmumps_ooc_init_vars dmumps_ooc_init_vars__
#define dmumps_ooc_end_write dmumps_ooc_end_write__
#elif defined(Add_)
#define dmumps_is_there_finished_request dmumps_is_there_finished_request_
#define dmumps_clean_request dmumps_clean_request_
#define dmumps_test_request dmumps_test_request_
#define dmumps_wait_request dmumps_wait_request_
#define dmumps_wait_all_requests dmumps_wait_all_requests_
#define dmumps_low_level_init_ooc_c dmumps_low_level_init_ooc_c_
#define dmumps_low_level_write_ooc_c dmumps_low_level_write_ooc_c_
#define dmumps_low_level_read_ooc_c dmumps_low_level_read_ooc_c_
#define dmumps_low_level_direct_read dmumps_low_level_direct_read_
#define dmumps_clean_io_data_c dmumps_clean_io_data_c_
#define dmumps_get_max_nb_req dmumps_get_max_nb_req_
#define dmumps_get_max_file_size dmumps_get_max_file_size_ 
#define dmumps_ooc_get_nb_files dmumps_ooc_get_nb_files_
#define dmumps_ooc_get_file_name dmumps_ooc_get_file_name_
#define dmumps_ooc_set_file_name dmumps_ooc_set_file_name_
#define dmumps_ooc_start_low_level dmumps_ooc_start_low_level_
#define dmumps_ooc_alloc_pointers dmumps_ooc_alloc_pointers_
#define dmumps_ooc_print_stats dmumps_ooc_print_stats_
#define dmumps_ooc_remove_file dmumps_ooc_remove_file_
#define dmumps_ooc_init_vars dmumps_ooc_init_vars_
#define dmumps_ooc_end_write dmumps_ooc_end_write_
#endif

#if defined(_WIN32)

//Next line May be needed depending on your Windows environment:
#define MUMPS_CALL __cdecl

#else
#define MUMPS_CALL
#endif

#define dmumps_ftnlen int

int MUMPS_CALL dmumps_is_there_finished_request(int* flag,int* ierr);

int MUMPS_CALL dmumps_clean_request(int* request_id,int* ierr);

int MUMPS_CALL dmumps_test_request(int* request_id,int* flag,int* ierr);

int MUMPS_CALL dmumps_wait_request(int* request_id,int* ierr);

int MUMPS_CALL dmumps_wait_all_requests(int* ierr);

int MUMPS_CALL dmumps_low_level_init_ooc_c(int* _myid, int* total_size_io,int* size_element,
                               int* async, int* k211, char* dmumps_dir, char* dmumps_file,
                               int* dmumps_dim_dir, int* dmumps_dim_file,
                               int* ierr, dmumps_ftnlen l1, dmumps_ftnlen l2);

int MUMPS_CALL dmumps_low_level_write_ooc_c( const int * strat_IO, 
                                 void * address_block,
                                 int * block_size,
                                 int * pos_in_file,
                                 int * file_number,
                                 int * inode,
                                 int * request_arg,
                                 int * ierr);

int MUMPS_CALL dmumps_low_level_read_ooc_c( const int * strat_IO, 
                                 void * address_block,
                                 int * block_size,
                                 int * from_where,
                                 int * file_number,
                                 int * inode,
                                 int * request_arg,
                                 int * ierr);

int MUMPS_CALL dmumps_low_level_direct_read(void * address_block,
                                 int * block_size,
                                 int * from_where,
                                 int * file_number,
                                 int * ierr);

int MUMPS_CALL dmumps_clean_io_data_c(int* myid,int* ierr);

int MUMPS_CALL dmumps_get_max_nb_req(int *max,int* ierr);

int MUMPS_CALL dmumps_get_max_file_size(double * max_ooc_file_size);

int MUMPS_CALL dmumps_ooc_get_nb_files(int* nb_files);

int MUMPS_CALL dmumps_ooc_get_file_name(int* indice,char* name,int* length, dmumps_ftnlen l1);

int MUMPS_CALL dmumps_ooc_set_file_name(int* indice,char* name,int* length,int* ierr, dmumps_ftnlen l1);

int MUMPS_CALL dmumps_ooc_alloc_pointers(int* dim,int* ierr);

int MUMPS_CALL dmumps_ooc_init_vars(int* myid_arg, int* nb_file_arg,
                                   int* size_element,int* async, int* k211,
                                   char* dmumps_dir, char* dmumps_file,
                                   int* dmumps_dim_dir, int* dmumps_dim_file,
                                   int *ierr, dmumps_ftnlen l1, dmumps_ftnlen l2);

int MUMPS_CALL dmumps_ooc_start_low_level(int* ierr);

int MUMPS_CALL dmumps_ooc_print_stats();

int MUMPS_CALL dmumps_ooc_remove_file(char *name,int* ierr, dmumps_ftnlen l1);

int MUMPS_CALL dmumps_ooc_end_write(int *ierr);

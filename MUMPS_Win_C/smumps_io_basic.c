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
/*    $Id: smumps_io_basic.c,v 1.79 2006/04/13 11:02:26 jylexcel Exp $  */

#include "smumps_io_basic_var.h"
#include "smumps_io_err_extern.h"

int _smumps_next_file(){
  /* Defines the pattern for the file name. The last 6 'X' will be replaced
     so as to name were unique */
  char name[150];
  int ret_code;


#ifndef _WIN32
  strcpy(name,smumps_ooc_file_prefix);
  mkstemp(name);
#else
  sprintf(name,"%s_%d",smumps_ooc_file_prefix,smumps_io_current_file_number+1);
#endif
  
  if (smumps_io_current_file_number >= smumps_io_nb_file-1){
    /* Exception : probably thrown because of a bad estimation
       of number of files. */
    /* We increase the number of file needed and then realloc. */
    smumps_io_nb_file++;
#ifndef _WIN32
    smumps_io_pfile_pointer_array=realloc((void *)smumps_io_pfile_pointer_array,smumps_io_nb_file*sizeof(int));
#else
    smumps_io_pfile_pointer_array=realloc((void *)smumps_io_pfile_pointer_array,smumps_io_nb_file*sizeof(FILE*));
#endif
    if(smumps_io_pfile_pointer_array==NULL){
      sprintf(error_str,"Allocation problem in low-level OOC layer\n");
      return -13;
    }
    smumps_io_pfile_name=realloc((void*)smumps_io_pfile_name,(smumps_io_nb_file)*sizeof(char *));
    if(smumps_io_pfile_name==NULL){
      sprintf(error_str,"Allocation problem in low-level OOC layer\n");
      return -13;
    }
  }
  


  /* Do change the current file */
  smumps_io_current_file_number++;
  /*  *(smumps_io_pfile_pointer_array+smumps_io_current_file_number)=fopen(name,"w+"); */
  *(smumps_io_pfile_name+smumps_io_current_file_number)=(char *)malloc((strlen(name)+1)*sizeof(char));

  if(*(smumps_io_pfile_name+smumps_io_current_file_number)==NULL){
    sprintf(error_str,"Allocation problem in low-level OOC layer\n");
    return -13;
  }
  strcpy(*(smumps_io_pfile_name+smumps_io_current_file_number),name);
  /* See smumps_io_basic.h for comments on the I/O flags passed to open */
#ifndef _WIN32
  *(smumps_io_pfile_pointer_array+smumps_io_current_file_number)=open(name,smumps_flag_open);
  if(*(smumps_io_pfile_pointer_array+smumps_io_current_file_number)==-1){
    smumps_io_build_err_str(errno,-90,"Unable to open OOC file",error_str,200);
    return -90;
  }
#else
  *(smumps_io_pfile_pointer_array+smumps_io_current_file_number)=fopen(name,"w");
  if(*(smumps_io_pfile_pointer_array+smumps_io_current_file_number)==NULL){
    /*    smumps_io_build_err_str(errno,-90,"Unable to open OOC file",error_str,200); */
    sprintf(error_str,"Problem while opening OOC file");
    return -90;
   }
#endif
  smumps_io_current_file=(smumps_io_pfile_pointer_array+smumps_io_current_file_number);
  
  smumps_io_current_file=(smumps_io_pfile_pointer_array+smumps_io_current_file_number);
  smumps_io_last_file_opened++;
  /*  if(*(smumps_io_pfile_pointer_array+smumps_io_current_file_number)==NULL) */
  smumps_io_write_pos=0;
  /*  printf("new file created %d\n",smumps_io_current_file_number);*/
  return 0;
}

void _smumps_update_current_file_position(){
    smumps_io_current_file_position=smumps_io_write_pos;
}


int _smumps_compute_where_to_write(const double to_be_written){
  /* Check if the current file has enough memory to receive the whole block*/
  double size;
  int ret_code;
  _smumps_update_current_file_position();
  size=(double)MAX_FILE_SIZE-((double)smumps_io_current_file_position+to_be_written);
  if (size > 0){
  }else{
    ret_code=_smumps_next_file();
    if(ret_code<0){
      return ret_code;
    }
    _smumps_update_current_file_position();
  }

  return 0;
}

int _smumps_prepare_pointers_for_write(double to_be_written,int * pos_in_file, int * file_number ){
  int ret_code;
  ret_code=_smumps_compute_where_to_write(to_be_written);
  if(ret_code<0){
    return ret_code;
  }
  *pos_in_file=smumps_io_current_file_position;
  *file_number=smumps_io_current_file_number;
  return 0;
}

int _smumps_compute_nb_concerned_files(int * block_size, int * nb_concerned_files){
  *nb_concerned_files=(int)ceil((((double)(*block_size))*((double)(smumps_elementary_data_size)))/((double)MAX_FILE_SIZE));
  return 0;
}

int smumps_io_do_write_block(void * address_block,
		     int * block_size,
		     int * pos_in_file,
		     int * file_number,		     
		     int * ierr){   
  /* Type of fwrite : size_t fwrite(const void *ptr, size_t size, 
                                    *size_t nmemb, FILE *stream); */
  size_t write_size;
  int i;
  int nb_concerned_files=0;
  int ret_code,file_number_loc,pos_in_file_loc;
  double to_be_written;
#ifndef _WIN32
  int* file;
#else
  FILE** file;
#endif
  int where;
  void* loc_addr;

  loc_addr=address_block;
  _smumps_compute_nb_concerned_files(block_size,&nb_concerned_files);
  to_be_written=((double)smumps_elementary_data_size)*((double)(*block_size));
  if((nb_concerned_files==0)&&(to_be_written==0)){
    *file_number=smumps_io_current_file_number;
    return 0;
  }

  for(i=0;i<nb_concerned_files;i++){

#ifndef _WIN32
#ifdef WITH_PFUNC
    if(smumps_io_flag_async==IO_ASYNC_TH){
      smumps_io_protect_pointers();
    }
#endif
#endif
    ret_code=_smumps_prepare_pointers_for_write(to_be_written,&pos_in_file_loc,&file_number_loc);
    
    if(ret_code<0){
#ifndef _WIN32
#ifdef WITH_PFUNC
    if(smumps_io_flag_async==IO_ASYNC_TH){
      smumps_io_unprotect_pointers();
    }
#endif
#endif
      return ret_code;
    }

    if(i==0){
      *pos_in_file=pos_in_file_loc;
      if(nb_concerned_files==1){
	*file_number=file_number_loc;
      }else{
	*file_number=-file_number_loc;
      }
    }
    if((double)(MAX_FILE_SIZE-smumps_io_write_pos)>to_be_written){
      write_size=(size_t)to_be_written;
    }else{
      write_size=(size_t)((double)(MAX_FILE_SIZE-smumps_io_write_pos));
    }
#if defined (_WIN32)
    write_size=(size_t)(int)((write_size)/smumps_elementary_data_size);
#endif
    file=smumps_io_current_file;
    where=smumps_io_write_pos;
#ifndef _WIN32
#ifdef WITH_PFUNC
    if(smumps_io_flag_async==IO_ASYNC_TH){
      smumps_io_unprotect_pointers();
    }
#endif
#endif
    /* printf("1 write -> size = %d | off = %d | file = %d\n",(int)write_size,where,*file_number);*/
    ret_code=smumps_io_write__(file,loc_addr,write_size,where);
    
    if(ret_code<0){
      return ret_code;
    }
    
#ifndef _WIN32
#ifdef WITH_PFUNC
    if(smumps_io_flag_async==IO_ASYNC_TH){
      smumps_io_protect_pointers();
    }
#endif
#endif

#ifndef _WIN32
    smumps_io_write_pos=smumps_io_write_pos+((int)write_size);
    to_be_written=to_be_written-((int)write_size);
    loc_addr=(void*)((size_t)loc_addr+write_size);
#else
    /* fread and write */
    smumps_io_write_pos=smumps_io_write_pos+((int)write_size*smumps_elementary_data_size);
    to_be_written=to_be_written-((int)write_size*smumps_elementary_data_size);
    loc_addr=(void*)((size_t)loc_addr+(size_t)((int)write_size*smumps_elementary_data_size));
#endif

#ifndef _WIN32
#ifdef WITH_PFUNC
    if(smumps_io_flag_async==IO_ASYNC_TH){
      smumps_io_unprotect_pointers();
    }
#endif
#endif
  }
  if(to_be_written!=0){
    sprintf(error_str,"Internal (1) error in low-level I/O operation %lf",to_be_written);
    return -90;
  }
  
  return 0;
}

int smumps_io_do_read_block(void * address_block,
	            int * block_size,
                    int * from_where,
		    int * file_number,
                    int * ierr){
  int ret_code;
  size_t size;
#ifndef _WIN32
  int* file;
#else
  FILE** file;
#endif
  double read_size;
  int local_fnum,local_offset;
  void *loc_addr;
  /*  if(((double)(*block_size))*((double)(smumps_elementary_data_size))>(double)MAX_FILE_SIZE){
    sprintf(error_str,"Internal error in low-level I/O operation (requested size too big for file system) \n");
    return -90;
    }*/

  if(*block_size==0)
    return 0;

  read_size=(double)smumps_elementary_data_size*(double)(*block_size);
  local_fnum=*file_number;
  if(local_fnum<0)
    local_fnum=-local_fnum;
  /*  if((*file_number<0)&&(read_size<(double)MAX_FILE_SIZE)){
    sprintf(error_str,"Internal error (1) in low level read op\n");
    return -90;
    }*/
  local_offset=*from_where;
  loc_addr=address_block;

  while(read_size>0){
    file=(smumps_io_pfile_pointer_array+(local_fnum));
#ifndef _WIN32
    if(read_size+(double)local_offset>(double)MAX_FILE_SIZE){
      size=(size_t)MAX_FILE_SIZE-(size_t)local_offset;
    }else{
      size=(size_t)read_size;
    }
#else
    if(read_size+(double)local_offset>(double)MAX_FILE_SIZE){
      size=((size_t)MAX_FILE_SIZE-(size_t)local_offset)/(size_t)smumps_elementary_data_size;
    }else{
      size=(size_t)(read_size/smumps_elementary_data_size);
    }
#endif
    ret_code=smumps_io_read__(file,loc_addr,size,local_offset);

    if(ret_code<0){
      return ret_code;
    }

#ifdef _WIN32
    size=size*smumps_elementary_data_size;
#endif

    read_size=read_size-(double)size;
    loc_addr=(void*)((size_t)loc_addr+size);
    local_fnum++;
    local_offset=0;
    if(local_fnum>smumps_io_nb_file){
      sprintf(error_str,"Internal error (2) in low level read op\n");
      return -90;
    }
  }
  return 0;
}
int smumps_free_file_pointers(){
  int i,ierr;
  for(i=0;i<=smumps_io_last_file_opened;i++){
#ifndef _WIN32
#ifdef SOLVE_CRAY
#endif
    ierr=close(*(smumps_io_pfile_pointer_array+i));
    if(ierr==-1){
      smumps_io_build_err_str(errno,-90,"Problem while closing OOC file",error_str,200);
      return -90;
    }
#else
    ierr=fclose(*(smumps_io_pfile_pointer_array+i));
    if(ierr==-1){
      sprintf(error_str,"Problem while opening OOC file\n");
      return -90;
    }    
#endif
    free(*(smumps_io_pfile_name+i));
  } 
  free(smumps_io_pfile_name);
  free(smumps_io_pfile_pointer_array);
  free(smumps_ooc_file_prefix);
  return 0;
}

int smumps_init_file_structure(int* _myid, int* total_size_io,int* size_element){
  /* Computes the number of files needed. Uses ceil value. */
  int ierr;
#ifndef _WIN32
  int k211_loc;
#endif

  smumps_io_nb_file=(int)((((double)(*total_size_io))*((double)(*size_element)))/(double)MAX_FILE_SIZE)+1;

  smumps_directio_flag=0;


  smumps_io_myid=*_myid;

  smumps_elementary_data_size=*size_element;
  /* Allocates the memory necessary to handle the file pointer array.*/
  ierr=smumps_io_alloc_file_struct(&smumps_io_nb_file);
  if(ierr<0){
    return ierr;
  }
  /* Init the current file.*/
  ierr=_smumps_next_file();
  if(ierr<0){
    return ierr;
  }
  return 0;
}
int smumps_init_file_name(char* smumps_dir,char* smumps_file,
			 int* smumps_dim_dir,int* smumps_dim_file,int* _myid){
  int i;
  char *tmp_dir,*tmp_fname;
  char base_name[20];
  int dir_flag=0,file_flag=0;
  char smumps_base[10]="smumps_";
  tmp_dir=(char *)malloc(((*smumps_dim_dir)+1)*sizeof(char));
  if(tmp_dir==NULL){
    sprintf(error_str,"Allocation problem in low-level OOC layer\n");
    return -13;
  }
  tmp_fname=(char *)malloc(((*smumps_dim_file)+1)*sizeof(char));
  if(tmp_fname==NULL){
    sprintf(error_str,"Allocation problem in low-level OOC layer\n");
    return -13;
  }
  for(i=0;i<*smumps_dim_dir;i++){
    tmp_dir[i]=smumps_dir[i];
  }
  tmp_dir[i]=0;
  for(i=0;i<*smumps_dim_file;i++){
    tmp_fname[i]=smumps_file[i];
  }
  tmp_fname[i]=0;  
  if(strcmp(tmp_dir,UNITIALIZED)==0){
    dir_flag=1;
    free(tmp_dir);
    tmp_dir=getenv("MUMPS_OOC_TMPDIR");
    if(tmp_dir==NULL){
#ifdef SP_
#ifndef CINES_
      tmp_dir=getenv("TMPDIR");
      if(tmp_dir==NULL){
	tmp_dir=SMUMPS_OOC_DEFAULT_DIR;
      }
#else
      tmp_dir=SMUMPS_OOC_DEFAULT_DIR;       
#endif       
#else
      tmp_dir=SMUMPS_OOC_DEFAULT_DIR;
#endif      
    }
  }
  if(strcmp(tmp_fname,UNITIALIZED)==0){
    free(tmp_fname);
    tmp_fname=getenv("MUMPS_OOC_PREFIX");
    file_flag=1;
  }

  if(tmp_fname!=NULL){
#ifndef _WIN32
      sprintf(base_name,"_%s%d_XXXXXX",smumps_base,*_myid);
#else
      sprintf(base_name,"_%s%d",smumps_base,*_myid);
#endif
      smumps_ooc_file_prefix=(char *)malloc((strlen(SEPARATOR)+strlen(tmp_dir)+strlen(tmp_fname)+strlen(base_name)+1+1)*sizeof(char));
      if(smumps_ooc_file_prefix==NULL){
	sprintf(error_str,"Allocation problem in low-level OOC layer\n");
	return -13;
      }
      sprintf(smumps_ooc_file_prefix,"%s%s%s%s",tmp_dir,SEPARATOR,tmp_fname,base_name);
  }else{
#ifndef _WIN32
    sprintf(base_name,"%s%s%d_XXXXXX",SEPARATOR,smumps_base,*_myid);
#else
    sprintf(base_name,"%s%s%d",SEPARATOR,smumps_base,*_myid);
#endif
      smumps_ooc_file_prefix=(char *)malloc((strlen(SEPARATOR)+strlen(tmp_dir)+strlen(base_name)+1)*sizeof(char));
      if(smumps_ooc_file_prefix==NULL){
	sprintf(error_str,"Allocation problem in low-level OOC layer\n");
	return -13;
      }
      sprintf(smumps_ooc_file_prefix,"%s%s%s",tmp_dir,SEPARATOR,base_name);
  }  
  if(!dir_flag){
    free(tmp_dir);
  }
  if(!file_flag){
    free(tmp_fname);
  }
  return 0;
}

int smumps_io_alloc_file_struct(int* nb){
#ifndef _WIN32  
  smumps_io_pfile_pointer_array=(int *)malloc((*nb)*sizeof(int));
  if(smumps_io_pfile_pointer_array==NULL){
    sprintf(error_str,"Allocation problem in low-level OOC layer\n");
    return -13;
  }
#else  
  smumps_io_pfile_pointer_array=(FILE **)malloc((*nb)*sizeof(FILE*));
  if(smumps_io_pfile_pointer_array==NULL){
    sprintf(error_str,"Allocation problem in low-level OOC layer\n");
    return -13;
  }
#endif

  smumps_io_pfile_name=(char **)malloc((*nb)*sizeof(char *));
  if(smumps_io_pfile_name==NULL){
#ifndef _WIN32      
    sprintf(error_str,"Allocation problem in low-level OOC layer\n");
#else
    sprintf(error_str,"Allocation problem in low-level OOC layer\n");
#endif    
    return -13;
    /*    return -1;*/
  }
  smumps_io_current_file_number = -1;
  smumps_io_last_file_opened=-1;
  return 0;
}

int smumps_io_get_nb_files(int* nb_files){
  *nb_files=smumps_io_last_file_opened+1;
  return 0;
}

int smumps_io_get_file_name(int* indice,char* name,int* length){
  int i;
  i=(*indice)-1;
  strcpy(name,*(smumps_io_pfile_name+i));
  *length=strlen(name)+1;
  return 0;  
}

int smumps_io_alloc_pointers(int * dim){
  return smumps_io_alloc_file_struct(dim);
}

int smumps_io_init_vars(int* myid_arg, int* nb_file_arg,int* size_element,int* async_arg){

#ifndef _WIN32
  int k211_loc;
#endif


  smumps_directio_flag=0;
  smumps_io_nb_file=*nb_file_arg;
  smumps_io_myid=*myid_arg;
  smumps_elementary_data_size=*size_element;
  smumps_io_flag_async=*async_arg;

  return 0;
}

int smumps_io_set_file_name(int* indice,char* name,int* length){
  int i;
  i=(*indice)-1;
  *(smumps_io_pfile_name+i)=(char *) malloc((*length)*strlen(name));
  if(*(smumps_io_pfile_name+i)==NULL){
    sprintf(error_str,"Allocation problem in low-level OOC layer");
    return -13;
  }
  strcpy(*(smumps_io_pfile_name+i),name);
  return 0;  
}

int smumps_io_open_files_for_read(){
  int i;

#ifdef IRIX64_
  struct dioattr dio;
#endif

  for(i=0;i<smumps_io_nb_file;i++){
#ifndef _WIN32
    *(smumps_io_pfile_pointer_array+i)=open(*(smumps_io_pfile_name+i),smumps_flag_open);
    
    if(*(smumps_io_pfile_pointer_array+i)==-1){
      smumps_io_build_err_str(errno,-90,"Problem while opening OOC file",error_str,200);
      return -90;
    }
#else
    *(smumps_io_pfile_pointer_array+i)=fopen(*(smumps_io_pfile_name+i),"a+");
    if(*(smumps_io_pfile_pointer_array+i)==NULL){
      sprintf(error_str,"Problem while opening OOC file");
      return -90;
    }
#endif
  }
  return 0;
}

int smumps_io_set_last_file(int* dim){
  smumps_io_last_file_opened=*dim-1;
  return 0;
}

#ifndef _WIN32  
#ifdef WITH_PFUNC

int smumps_io_protect_pointers(){
  pthread_mutex_lock(&smumps_io_pwrite_mutex);
  return 0;
}
int smumps_io_unprotect_pointers(){
  pthread_mutex_unlock(&smumps_io_pwrite_mutex);
  return 0;
}

int smumps_io_init_pointers_lock(){
  pthread_mutex_init(&smumps_io_pwrite_mutex,NULL);
  return 0;
}
int smumps_io_destroy_pointers_lock(){
  pthread_mutex_destroy(&smumps_io_pwrite_mutex);
  return 0;
}

#endif /*WITH_PFUNC*/
#endif /*_WIN32*/

int smumps_io_read__(void * file,void * loc_addr,size_t size,int local_offset){
  int ret_code;
#ifndef _WIN32
  if(!smumps_directio_flag){
    ret_code=smumps_io_read_os_buff__(file,loc_addr, size,local_offset);
    if(ret_code<0){
      return ret_code;
    }
  }
#else
  ret_code=smumps_io_read_win32__(file,loc_addr, size,local_offset);
  if(ret_code<0){
    return ret_code;
  }
#endif  
  return 0;
}

#ifndef _WIN32
int smumps_io_read_os_buff__(void * file,void * loc_addr,size_t size,int local_offset){
  size_t ret_code;
  /*  printf("Read with buff %d %d %d\n",(int) size, local_offset,*((int *)file));*/
#ifdef WITH_PFUNC
  ret_code=pread(*(int *)file,loc_addr,size,local_offset);
#else
  lseek(*(int *)file,(long) local_offset,SEEK_SET);
  ret_code=read(*(int *)file,loc_addr,size);
#endif
  if((int) ret_code==-1){
    smumps_io_build_err_str(errno,-90,"Problem with low level read",error_str,200);
    return -90;
  }
  return 0;
}

#endif

#ifdef _WIN32
int smumps_io_read_win32__(void * file,void * loc_addr,size_t size,int local_offset){
  size_t ret_code;
  fseek(*(FILE **)file,(long) local_offset,SEEK_SET);
  ret_code=fread(loc_addr,smumps_elementary_data_size,size,*(FILE **)file);
  if(ret_code!=size){
    sprintf(error_str,"Problem with I/O operation\n");
    return -90;
  }
  return 0;
}
#endif

int smumps_io_write__(void *file, void *loc_addr, size_t write_size, int where){
  int ret_code;
#ifndef _WIN32
  if(!smumps_directio_flag){
    ret_code=smumps_io_write_os_buff__(file,loc_addr, write_size,where);
    if(ret_code<0){
      return ret_code;
    }
  }
#else
  ret_code=smumps_io_write_win32__(file,loc_addr, write_size,where);
  if(ret_code<0){
    return ret_code;
  }
#endif
  return 0;
}

#ifndef _WIN32
int smumps_io_write_os_buff__(void *file, void *loc_addr, size_t write_size, int where){
  size_t ret_code;
#ifdef WITH_PFUNC
  ret_code=pwrite(*(int *)file,loc_addr,write_size,where);
#else
  /*in this case all the I/O's are made by the I/O thread => we don't
    need to protect the file pointer.*/
  lseek(*(int *)file,(long)where,SEEK_SET); 
  ret_code=write(*(int *)file,loc_addr,write_size);
#endif
  if(((int)ret_code==-1)||(ret_code!=write_size)){
    smumps_io_build_err_str(errno,-90,"Problem with low level write",error_str,200);
    return -90;
  }
  return 0;
}

#endif

#ifdef _WIN32
int smumps_io_write_win32__(void *file, void *loc_addr, size_t write_size, int where){
  size_t ret_code;
  fseek(*(FILE **)file,(long)where,SEEK_SET);  
  ret_code=fwrite(loc_addr,smumps_elementary_data_size, write_size,*(FILE**)file);
  if(ret_code!=write_size){
    sprintf(error_str,"Problem with I/O operation\n");
    return -90;
  }
  return 0;
}
#endif



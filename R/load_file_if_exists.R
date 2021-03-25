load_file_if_exists = function(path){
 if(file.exists(path)){
 
   load(path, envir=.GlobalEnv)
 }
  
}


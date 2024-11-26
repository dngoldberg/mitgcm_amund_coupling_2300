function write_coupled_input(fname,arr,fold)

global sub_dir

 fname2 = [sub_dir '/' fold '/' fname];
 binwrite(fname2,arr);

end

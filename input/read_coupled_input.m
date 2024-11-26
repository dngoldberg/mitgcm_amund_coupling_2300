function q = read_coupled_input(fname,fold,varargin)

global sub_dir

 fname2 = [sub_dir '/' fold '/' fname];
 if(length(varargin)==1);
  q = binread(fname2,8,varargin{1});
 elseif(length(varargin)==2);
  q = binread(fname2,8,varargin{1},varargin{2});
 elseif(length(varargin)==3);
  q = binread(fname2,8,varargin{1},varargin{2},varargin{3});
 end 

end

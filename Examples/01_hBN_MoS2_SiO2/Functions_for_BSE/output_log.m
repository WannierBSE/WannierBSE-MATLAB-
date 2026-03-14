function output_log(DataPath,N_kpt,Nv,Nc,Npar_HBSE,Time)
        

cd(DataPath)
%----------------------------------------------------------
% record in computational conditions
fid = fopen('WBSE_log.txt','w');

t_now = datetime;
fprintf(fid,'This calculation was finished at %s \n \n',t_now);

fprintf(fid,'-------------- k-points and bands and included in BSE --------------  \n');
fprintf(fid,'number of k-points Nk = %d \n',N_kpt);
fprintf(fid,'numbers of valence and conduction bands (Nv,Nc) = ( %d , %d ) \n \n \n',Nv,Nc);

fprintf(fid,'-------------- CPU usage -------------- \n');
fprintf(fid,'The number of cores used in the construction of BSE matrix = %d \n \n',Npar_HBSE);

fprintf(fid,'-------------- Computation time -------------- \n');
fprintf(fid,'Computation time of constructing  BSE matrix  =   %12.3f sec. \n',Time(1));
fprintf(fid,'Computation time of diagonalizing BSE matrix  =   %12.3f sec. \n',Time(2));
fprintf(fid,'                      Total computation time  =   %12.3f sec. \n',Time(3));


fclose(fid);

end
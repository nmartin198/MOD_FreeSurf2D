function rIToUTPUT
    % This script is a component of "MOD_FreeSurf2D: a Matlab surface
    %   fluid flow model for rivers and streams."
    %   by N. Martin and S. Gorelick (2004)
    %
    % rIToUTPUT writes the fluid flow information out to matlab .mat files
    % at to .txt files.

    global H Hux Hvy u v XINDEX YINDEX

    % local variables.
    % output.
    File1 = fopen('U.txt','w+');
    File2 = fopen('V.txt','w+');
    File3 = fopen('Hux.txt','w+');
    File4 = fopen('Hvy.txt','w+');
    File5 = fopen('H.txt','w+');
    File6 = fopen('XINDEX.txt','w+');
    File7 = fopen('YINDEX.txt','w+');
    fprintf(File1,'%10.6f\n',u);
    fprintf(File2,'%10.6f\n',v);
    fprintf(File3,'%10.6f\n',Hux);
    fprintf(File4,'%10.6f\n',Hvy);
    fprintf(File5,'%10.6f\n',H);
    fprintf(File6,'%10.6f\n',XINDEX);
    fprintf(File7,'%10.6f\n',YINDEX);
    fclose(File1);
    fclose(File2);
    fclose(File3);
    fclose(File4);
    fclose(File5);
    fclose(File6);
    fclose(File7);
    % not needed
    %clear File1 File2 File3 File4 File5 File6 File7;
end
%EOF
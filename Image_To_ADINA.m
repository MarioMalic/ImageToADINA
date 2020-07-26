%%
% Mario Malic, July 2020
% Contact mail: mario@mariomalic.com
% Alternative contact mail: mario.malic@gmail.com
%
% Credits:
% - Function 'rgb2hex' written by Chad A. Greene, sourced from
%   MATLAB File Exchange
%% Recommendations:
% There are few options to reduce the time building the model
% - Edit - Memory Usage and give ADINA-AUI extra RAM to work with.
% - File - and click Stop Recording if it is turned on
% - Hide message window on the bottom of AUI window
% On my machine 161 x 240 image takes 10 minutes, 402 x 600 takes an hour
% 515 x 768 six hours
% ADINA will most likely freeze (not respond) during the model buildup, have patience.
%% Issues:
% - RGB colors are not well represented.
%% Possible solutions:
% - Check if RGB_Array, Deformation_Array have been correctly assigned/calculated
% - Investigate Band_Plot capabilities
%% Postprocessing the model
% - This script outputs Bandtable_Unique.txt file that contains the values
%   custom for table of strain values and corresponding RGB values in HEX format.
% - If you use the Type: Automatic band plot, Min/Max color is the first
%   and last HEX value in the mentioned text file or from variables R_Min,
%   G_Min, B_Min, R_Max, G_Max, B_Max 
% - One can use black and white as Min/Max values for a grayscale image, looks great
%%
clearvars;
clc;
fclose all;
format long
%% Obtaining data for the image
% Image path (Next 2 lines need adjustment)
Working_Dir = 'C:\Users\Mario\Desktop\ADINA\Last_Test\'; % Folder where the image is located in, including the last character '\'
File_Name = 'Mona_Lisa_Smallest'; % Image name
cd (Working_Dir);
Image_Data = imread (File_Name, 'jpg'); % Works only with .jpg format

% RGB Channels
Red_Channel = single(Image_Data(:, :, 1));
Green_Channel = single(Image_Data(:, :, 2));
Blue_Channel = single(Image_Data(:, :, 3));

% Horizontal and vertical pixels
[V_Pixels, H_Pixels, ~] = size(Image_Data);

% Contains pixel data as a single value
RGB_Array =double(65536 * Red_Channel + 256 * Green_Channel + Blue_Channel);

% Finding color with lowest and highest RGB value
Min_RGB = min(min(RGB_Array));
[X_Min, Y_Min] = find(RGB_Array == Min_RGB);
Max_RGB = max(max(RGB_Array));
[X_Max, Y_Max] = find(RGB_Array == Max_RGB);

% Finding minimum and maximum RGB values
% These values are the input for the minimum and maximum color for the band plot
R_Min = Red_Channel(X_Min, Y_Min)
G_Min = Green_Channel(X_Min, Y_Min)
B_Min = Blue_Channel(X_Min, Y_Min)
R_Max = Red_Channel(X_Max, Y_Max)
G_Max = Green_Channel(X_Max, Y_Max)
B_Max = Blue_Channel(X_Max, Y_Max)

% Calculate needed deformations for each pixel
Deformation_Array = rescale(RGB_Array, 1e-6, 1);

% Obtaining values for custom band table
Values = cell(V_Pixels, H_Pixels);
Colors = cell(V_Pixels, H_Pixels);
Counter = 1;
Band_Table = cell(V_Pixels*H_Pixels, 2);

for ii = 1 : 1 : V_Pixels
    for jj = 1 : 1 : H_Pixels 
        Values {ii, jj} = Deformation_Array(ii,jj);
        Colors {ii, jj} = rgb2hex([Red_Channel(ii,jj) Green_Channel(ii,jj) Blue_Channel(ii,jj)]);
        Band_Table {Counter, 1} = sprintf('%.16f\t', Values {ii, jj});
        Band_Table {Counter, 2} = sprintf('%s', Colors {ii, jj});
        Counter = Counter + 1;
    end
end

Band_Table_New = Band_Table';
Band_Table_1 = Band_Table_New(1,:);
Band_Table_2 = Band_Table_New(2,:);
Band_Table_New= strcat(Band_Table_1, Band_Table_2);

Band_Table = (sort(Band_Table))';
Band_Table_Unique_1 = unique(Band_Table(1,:));
Band_Table_Unique_2 = unique(Band_Table(2,:));
Band_Table_Unique = strcat(Band_Table_Unique_1, Band_Table_Unique_2);

% Writing custom band table
fid = fopen('Bandtable.txt','w');
fprintf(fid,'%s\n', Band_Table_New{:});
fclose(fid);

% Writing custom band table for direct import
fid = fopen('Bandtable_Unique.txt','w');
fprintf(fid,'%s\n', Band_Table_Unique{:});
fclose(fid);
%% Writing input file for ADINA

DX = 1e-4; % Distance between horizontal pixels
DY = 1e-4; % Distance between vertical pixels
Num_Points = H_Pixels * V_Pixels; 
CS_Array = zeros(V_Pixels, H_Pixels); % Coordinate system for points input

X_Coords_Array = zeros(V_Pixels, H_Pixels);
for ii = 1 : 1 : V_Pixels
    for jj = 1 : 1 : H_Pixels   
      X_Coords_Array(ii,jj) = DX*(jj-1); 
    end
end

Y_Coords_Array = zeros(V_Pixels, H_Pixels);
for ii = 1 : 1 : V_Pixels
    for jj = 1 : 1 : H_Pixels   
      Y_Coords_Array(ii,jj) = (-1)*DY*(ii-1); % Function "imread" reads pixels in a different orientation, hence the first (-1)
    end
end

Z_Coords_Array_0 = ones(V_Pixels, H_Pixels)*(-1);
Z_Coords_Array_1 = ones(V_Pixels, H_Pixels) - Deformation_Array;

%% Writing to file

% Changing environment settings to speed up the model build up
Index = 1;
File{Index} = sprintf('CONTROL PLOTUNIT=PERCENT VERBOSE=YES ERRORLIM=0 LOGLIMIT=0 UNDO=-1 PROMPTDE=UNKNOWN AUTOREPA=NO DRAWMATT=NO DRAWTEXT=EXACT DRAWLINE=EXACT DRAWFILL=EXACT AUTOMREB=NO ZONECOPY=NO SWEEPCOI=YES SESSIONS=NO DYNAMICT=YES UPDATETH=YES AUTOREGE=NO ERRORACT=CONTINUE FILEVERS=CURRENT INITFCHE=NO SIGDIGIT=12 AUTOZONE=NO PSFILEVE=V0 ELEMENT-=REPEAT SESSIONO=ALL SESSIONT=100 SZNAME=TOPDOWN DUPLICAT=NO');
Index = Index + 1; File{Index} = sprintf('*');

% Coordinate system header
Index = Index + 1;File{Index} = sprintf('COORDINATES POINT SYSTEM=0');
Index = Index + 1; File{Index} = sprintf('@CLEAR');

% Defining points for set 0 (Z = -1)
Point_Index = 0;
for ii = 1 : 1 : V_Pixels
    for jj = 1 : 1 : H_Pixels
        Index = Index + 1; % Starts from 1 in loop
        Point_Index = Point_Index + 1;
        File{Index} = sprintf('%d %.16f %.16f %.16f %d', Point_Index, X_Coords_Array(ii,jj), Y_Coords_Array(ii,jj), Z_Coords_Array_0(ii,jj), CS_Array(ii,jj));
    end
end

% Defining points for set 1 (Z = 1 - Deformation_Array)
for ii = 1 : 1 : V_Pixels
    for jj = 1 : 1 : H_Pixels
        Index = Index + 1; 
        Point_Index = Point_Index + 1;
        File{Index} = sprintf('%d %.16f %.16f %.16f %d', Point_Index, X_Coords_Array(ii,jj), Y_Coords_Array(ii,jj), Z_Coords_Array_1(ii,jj), CS_Array(ii,jj));
    end
end
Index = Index + 1; File{Index} = sprintf('@');
Index = Index + 1; File{Index} = sprintf('*');

% Defining lines
Line_Index = 0;
for ii = 1 : Num_Points
    Line_Index = Line_Index + 1;
    Index = Index + 1;
    File{Index} = sprintf('LINE STRAIGHT NAME=%d P1=%d P2=%d', Line_Index, ii, Num_Points + ii);
end
Index = Index + 1; File{Index} = sprintf('*');

% Defining boundary conditions
Index = Index + 1; File{Index} = sprintf('FIXBOUNDARY POINTS FIXITY=ALL'); % BC Header
Index = Index + 1; File{Index} = sprintf('@CLEAR'); 
for ii = 1 : Num_Points
    Index = Index + 1;
    File{Index} = sprintf('%d ''ALL''', ii);
end
Index = Index + 1; File{Index} = sprintf('@');
Index = Index + 1; File{Index} = sprintf('*');


% Defining load (displacement)
Load_Array = Deformation_Array;

Load_Index = 0;
for ii = 1 : 1 : V_Pixels
    for jj = 1 : 1 : H_Pixels
        Load_Index = Load_Index + 1;
        Index = Index + 1;
        File{Index} = sprintf('LOAD DISPLACEMENT NAME=%d DX=FREE DY=FREE DZ=%.16f AX=FREE, AY=FREE AZ=FREE', Load_Index, Load_Array(ii,jj));
    end
end
Index = Index + 1; File{Index} = sprintf('*');

% Applying load (displacement)
Index = Index + 1; File{Index} = sprintf('APPLY-LOAD BODY=0');
Index = Index + 1; File{Index} = sprintf('@CLEAR');

Load_Index = 0;
for ii = Num_Points + 1 : 1 : 2*Num_Points
    Load_Index = Load_Index + 1;
    Index = Index + 1;
    File{Index} =  sprintf('%d  ''DISPLACEMENT'' %d  ''POINT'' %d 0 1 0.00000000000000 0 -1 0 0 0 ''NO'' 0.00000000000000 0.00000000000000 1 0  ''MID''', Load_Index, Load_Index, ii);
end
Index = Index + 1; File{Index} = sprintf('@');
Index = Index + 1; File{Index} = sprintf('*');

% Defining material model
Index = Index + 1; File{Index} = sprintf('MATERIAL ELASTIC NAME=1 E=2.10000000000000E+11 NU=0.300000000000000 DENSITY=0.00000000000000 ALPHA=0.00000000000000 MDESCRIP=''NONE''');
Index = Index + 1; File{Index} = sprintf('*');

% % Defining cross section
% Index = Index + 1; File{Index} = sprintf('CROSS-SECTIO RECTANGULAR NAME=1 WIDTH=0.000100000000000000 HEIGHT=0.000100000000000000 SC=0.00000000000000 TC=0.00000000000000 TORFAC=1.00000000000000 SSHEARF=0.00000000000000 TSHEARF=0.00000000000000 ISHEAR=NO SQUARE=NO');
% Index = Index + 1; File{Index} = sprintf('*');

% Defining element group
Index = Index + 1; File{Index} = sprintf('EGROUP TRUSS NAME=1 SUBTYPE=GENERAL DISPLACE=DEFAULT MATERIAL=1 INT=DEFAULT GAPS=NO INITIALS=NONE CMASS=DEFAULT TIME-OFF=0.00000000000000 OPTION=NONE RB-LINE=1 DESCRIPT=NONE AREA=1.00000000000000 PRINT=DEFAULT SAVE=DEFAULT TBIRTH=0.00000000000000 TDEATH=0.00000000000000 TMC-MATE=1 RUPTURE=ADINA GAPWIDTH=0.00000000000000');
Index = Index + 1; File{Index} = sprintf('*');

% Meshing lines
Index = Index + 1; File{Index} = sprintf('GLINE NODES=2 NCOINCID=NO SUBSTRUC=0 GROUP=1 MIDNODES=CURVED XO=0.00000000000000 YO=0.00000000000000 ZO=0.00000000000000 XYZOSYST=SKEW');
Index = Index + 1; File{Index} = sprintf('@CLEAR');
for ii = 1 : 1 : Num_Points
    Index = Index + 1;
    File{Index} = sprintf('%d', ii);
end
Index = Index + 1; File{Index} = sprintf('@');
Index = Index + 1; File{Index} = sprintf('*');

% Defining DOFs
Index = Index + 1; File{Index} = sprintf('MASTER ANALYSIS=STATIC MODEX=EXECUTE TSTART=0.00000000000000 IDOF=110111 OVALIZAT=NONE FLUIDPOT=AUTOMATIC CYCLICPA=1 IPOSIT=STOP REACTION=YES INITIALS=NO FSINTERA=NO IRINT=DEFAULT CMASS=NO SHELLNDO=AUTOMATIC AUTOMATI=ATS SOLVER=SPARSE CONTACT-=CONSTRAINT-FUNCTION TRELEASE=0.0000000000000 RESTART-=NO FRACTURE=NO LOAD-CAS=NO LOAD-PEN=NO SINGULAR=YE STIFFNES=0.000100000000000000 MAP-OUTP=NONE MAP-FORM=N NODAL-DE='''' POROUS-C=NO ADAPTIVE=0 ZOOM-LAB=1 AXIS-CYC= PERIODIC=NO VECTOR-S=GEOMETRY EPSI-FIR=NO STABILIZ=NO STABFACT=1.00000000000000E-10 RESULTS=PORTHOLE FEFCORR=NO BOLTSTEP=1 EXTEND-S=YES CONVERT-=NO DEGEN=YES TMC-MODE=NO ENSIGHT-=NO IRSTEPS=1 INITIALT=NO TEMP-INT=NO ESINTERA=NO OP2GEOM=NO INSITU-D=NO OP2ERCS=ELEMENT 2DPL-AX=YZ-Z OP2STR=DEFAULT IRTIMES=0');
Index = Index + 1; File{Index} = sprintf('*');

% Defining kinematic assumptions
Index = Index + 1; File{Index} = sprintf('KINEMATICS DISPLACE=LARGE STRAINS=SMALL UL-FORMU=DEFAULT PRESSURE=NO INCOMPAT=AUTOMATIC RIGIDLIN=NO BEAM-ALG=CURRENT KBEAM-EI=NO BACKSTRE=NO');
Index = Index + 1; File{Index} = sprintf('*');

% Saving file in .IDB format (not recommended due to big file size)
% Index = Index + 1; File{Index} = sprintf('DATABASE SAVE PERMFILE=''%s%s.idb'' PROMPT=NO', Working_Dir, File_Name);
% Index = Index + 1; File{Index} = sprintf('*');

% Exporting file in .NAS format 
Index = Index + 1; File{Index} = sprintf('EXPORT NASTRAN FILE=''%s%s.nas'' OVERWRIT=YES FORMAT=SMALL', Working_Dir, File_Name);
Index = Index + 1; File{Index} = sprintf('*');

% Saving input file
Input_File = sprintf('%s.in', File_Name);
FID = fopen(Input_File,'wt'); % Saves .in file
for p = 1:numel(File)
    fprintf(FID,'%s\n', File{p});
end
fclose(FID);


function [ flag ] = bootstrap_hypoDD( path,path2,boot_num,verbose,HYPODD)
%% Arguments description:
%INPUT  path: where the original hypoDD relocation files are located
%       main_path: bootstrap working directory.
%       boot_num: the number of 
%       verbose: 1, run time display on; 0. run time display off.
%       HYPODD: path of command HYPODD.
%OUTPUT flag: flag = 1, code runs successfully;  flag = 0, failed.
%% Setting parameters
new_dt = 'new_dt'; %For dt.ct files created randomly by this function
new_reloc = 'new_reloc';% For newly created relocation files 
out = 'bootstrap.out';% Output file containing bootstrap errors.
c_path = pwd;
verbose = 1;
data = [];
flag = 0;
% MAXDATA = 300000;% Max number of dtimes;
%% Check if directories exist 
if verbose == 1
disp('Check if directories exist...');
end
if exist(path2) ~=7
    mkdir(path2);
end
cd(path2)
if exist(path) ~=7
    error(['Relocation directory ','''','path','''',' doesn not exist!']);
elseif exist([path,'/dt.ct']) == 0 | exist([path,'/hypoDD.res']) == 0
    error('''dt.ct'' or ''hypoDD.res'' does not exist.');
end 

%% Creating working directories in function running folder
if verbose == 1
disp('Creating working forectories...');
end

if exist(['./',new_dt])  == 7
    rmdir(['./',new_dt],'s')
end
if exist(['./',new_reloc]) == 7
    rmdir(['./',new_reloc],'s')
end
    mkdir(new_dt)
    mkdir(new_reloc)
%% get TT info; initial file
if verbose == 1
disp('Reading dt.ct...');
end
fid_dt = fopen([path,'/dt.ct']);
c = 0;
j = 0;
num_p = 0;
num_s = 0;
while ~feof(fid_dt)
    line2 = fgetl(fid_dt);
    j = j + 1;
    if (line2 < 0)
        break;
    end
    
    if (sscanf(line2(1), '%s') == '#') %go to next line if this is the case
        c = c+1;
        C1 = sscanf(line2(8:11), '%d');
        C2 = sscanf(line2(18:21), '%d');
        line2 = fgetl(fid_dt);
    end
    
    iC1(j) = C1; iC2(j) = C2;
    sta(j) = {sscanf(line2(1:4), '%s')};
    TT1(j) = sscanf(line2(10:15), '%f');
    TT2(j) = sscanf(line2(18:23), '%f');
    wt(j) = sscanf(line2(25:30), '%f');
    phase(j) = {sscanf(line2(length(line2)),'%s')};
    if char(phase(j)) == 'P'
        num_p = num_p+1;
    elseif char(phase(j)) == 'S'
        num_s = num_s + 1;
    end
end
fclose(fid_dt);p_type = char(phase);
S_phase_index = find(p_type=='S');
P_phase_index = find(p_type=='P');


%% get hypoDD.res
if verbose == 1
disp('Reading hyoDD.res...');
end
H_res = dlmread([path,'/hypoDD.res'],'',1,1);
ID_res = H_res(:,4);
res = H_res(:,6);
res_P = res(find(H_res(:,4)==3));
res_S = res(find(H_res(:,4)==4));
    if length(res_S) <= 1
     res_S = zeros(1);
     disp('# of res_S <= 1!');
    end
    if length(res_P) <= 1
     res_P = zeros(1);
     disp('# of res_P <= 1!');

    end

%% randomly add residual information to TT1
if verbose == 1
disp('Generating dtrand files...')
end
for c = 1 :boot_num
    vecRes_P = randperm(length(res_P)); %random indeces for res_P
    vecRes_S = randperm(length(res_S)); % random indecs for res_S
    vecPindex = randperm(length(P_phase_index));
    vecSindex = randperm(length(S_phase_index));
    vecTTP = P_phase_index(vecPindex); %random indeces for TTP
    vecTTP = vecTTP(1:length(vecRes_P)); %truncate random indeces to length of residual vector
    vecTTS = S_phase_index(vecSindex); %random indeces for TTS
    vecTTS = vecTTS(1:length(vecRes_S)); %truncate random indeces to length of residual vector

    perturbT = TT1; %so that this new file can be perturbed, and not sacrifice the old one
    for k = 1:length(res_P)
        indTT = vecTTP(k);
        perturbT(indTT) = TT1(vecTTP(k)) + (res_P(vecRes_P(k))/1000);
    end
    for k = 1:length(res_S)
        indTT = vecTTS(k);
        perturbT(indTT) = TT1(vecTTS(k)) + (res_S(vecRes_S(k))/1000);
        res_S(vecRes_S(k));
    end

 %now to put back into a dt.ct file to save for analysis
   lb = '#';
    dtfile = strcat('dtrand',num2str(c),'.ct');
    if verbose == 1
    fprintf('Writing %s\n',dtfile);
    end
    fid = fopen([new_dt,'/',dtfile], 'wt');

    for q = 1:length(TT1)
        if (q == 1)
            fprintf(fid, '%s       %3d       %3d\n', lb, iC1(q), iC2(q));
            fprintf(fid, '%s     %6.3f  %6.3f %6.4f %s\n', char(sta(q)), perturbT(q), TT2(q), wt(q), char(phase(q)));
        elseif (q > 1 && iC1(q) == iC1(q-1) && iC2(q) == iC2(q-1))
            fprintf(fid, '%s     %6.3f  %6.3f %6.4f %s\n', char(sta(q)), perturbT(q), TT2(q), wt(q), char(phase(q)));
        elseif (q > 1 && (iC1(q) ~= iC1(q-1) || iC2(q) ~= iC2(q-1)))
            fprintf(fid, '%s       %3d       %3d\n', lb, iC1(q), iC2(q));
            fprintf(fid, '%s     %6.3f  %6.3f %6.4f %s\n', char(sta(q)), perturbT(q), TT2(q), wt(q), char(phase(q)));
        end
    end
    fclose(fid);
    tmp = sprintf('%03d',c);
    if exist(['./',new_reloc]) ~=7
    mkdir([new_reloc,'/',tmp]);
    end
    copyfile([path,'/*'],[new_reloc,'/',tmp]);
    copyfile([new_dt,'/',dtfile],[new_reloc,'/',tmp,'/dt.ct']);
end

%% Relocate using new dt files 
for c = 1:boot_num
    tmp = sprintf('%03d',c);
    cd([new_reloc,'/',tmp]);
    system([HYPODD,' ./hypoDD.inp']);
    M = load('hypoDD.reloc');
    data(end+1:end+size(M,1),:) = M(:,1:4);
    cd(c_path);
end
ID_m = data(:,1); ex=data(:,3);ey=data(:,2);ez=data(:,4);
%% Read original relocation file
m = load([path,'hypoDD.reloc']); id_o = m(:,1);

%% Evaluate boostrap errors
fid= fopen('error_analysis.out','w');
fprintf(fid,'   LAT         LON        DEP     AX1(km) AX2(km) AY1(km) AY2(km) AZ1(km) AZ2(km) EX(km) EY(km)  EZ(km)  YR   MO DY HR MN  SC    MG    ID\n');
fprintf('   LAT         LON        DEP     AX1(km) AX2(km) EZ(km) AZ  YR   MO DY HR MN  SC    MG    ID\n');
for i=1:length(id_o)
        ii= find(ID_m==id_o(i));
	if(length(ii)>10) % && id(i) ~= 60);
        [x,y,z] = lla2xyz(mean(ey(ii)),mean(ex(ii)),ey(ii),ex(ii),ez(ii));
        [ ax1,ax2]= error_ellipse_fnc(x,y,0);
        [ay1,ay2]= error_ellipse_fnc(x,z,0);
        [az1,az2] = error_ellipse_fnc(y,z,0);
        stdx=std(x);stdy=std(y);stdz=std(z);
        fprintf(fid,'%11.6f%12.6f%8.3f%8.3f%8.3f%8.3f%8.3f%8.3f%8.3f%8.3f%8.3f%8.3f%5i%3i%3i%3i%3i%7.3f%5.1f%10i\n',m(i,2),m(i,3),m(i,4),ax1,ax2,ay1,ay2,az1,az2,...
                               stdx,stdy,stdz,m(i,11),m(i,12),m(i,13),m(i,14),m(i,15),m(i,16),m(i,17),m(i,1));
  else
        ax1=-9; ax2=-9;
        e2=-9; az=-9;
	end
end; 
fclose(fid);

%% Return flag
if exist([path2,'/bootstrap_done'])
   delete([path2,'/bootstrap_done']) ;
end
system(['touch ',path2,'/bootstrap_done!']);
flag = 1; 
cd(c_path);
return

end


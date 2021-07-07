% Run several METiS 
clear all
close all
clc

if ispc
    folder = 'G:\Meu Drive\Doutorado\1Testes_Jappaku\metis\BIC\180DEG\CDz_real';
    
    % Remove files from previous runs. This is only needed because the next step
    % is to list all the .txt files in the folder.
    system(['del "' folder '\*_log*.txt"']);
    system(['del "' folder '\*_sum*.txt"']);
    system(['del "' folder '\*_out*.txt"']);

    files_to_run = dir([folder '\*.txt']);
    
    % Run each of the files in sequence
    for ii = 1 : numel(files_to_run)
        tic
        flNm = [files_to_run(ii).folder '\' files_to_run(ii).name];
        [status, ~] = system(['"C:\Users\lucas.henrique\Documents\METiS - VS\x64\Release\METiS - VS.exe" "' flNm '"']);
%         disp(['status:' status]);
        fprintf('%d files left - Last file took %.3f s\n', numel(files_to_run) - ii, toc);
    end



% % % % % elseif isunix
% % % % %     folder = './metis_OC4_completo_h200_45deg_CD0/A_1m';
% % % % %     
% % % % %     % Remove files from previous runs. This is only needed because the next step
% % % % %     % is to list all the .txt files in the folder.
% % % % %     system(['rm ' folder '/*_log*.txt']);
% % % % %     system(['rm ' folder '/*_sum*.txt']);
% % % % %     system(['rm ' folder '/*_out*.txt']);
% % % % %     
% % % % %     files_to_run = dir([folder '/*.txt']);
% % % % %     
% % % % %     % Run each of the files in sequence
% % % % % %     for ii = 1 : numel(files_to_run)
% % % % % %         flNm = [files_to_run(ii).folder '\' files_to_run(ii).name];
% % % % % %         system(['"Z:\VS Project - METiS\x64\Release\METiS.exe" "' flNm '"']);
% % % % % %     end
% % % % %     
% % % % %    

end


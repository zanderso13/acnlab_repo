datadir = '/Volumes/ZachExternal/ACNlab/MWMH/dti/';
session = 'ses-1';
tracks = {'Uncinate','LAmyg_mOFC','LNAcc_Amyg','LNAcc_mOFC','RAmyg_mOFC','RNAcc_Amyg','RNAcc_mOFC'};

%% Snag those files
for trk = 1:length(tracks)
    
    % extract specific stats after opening the files. The format of the
    % files is not super conducive to this software. Also, the place where
    % numbers appear differs depending on what options you use. 
    if strcmp(tracks{trk}, 'Uncinate') == 1
        fnames_left = filenames(fullfile(datadir,session,'results',tracks{trk},'*left*.stat.txt'));
        fnames_right = filenames(fullfile(datadir,session,'results',tracks{trk},'*right*.stat.txt'));
        for sub = 1:length(fnames_left)   
            fid = fopen(fnames_left{sub});
            datal = textscan(fid,'%s%s%s');
            fclose(fid);

            num_tracts_left(sub,1)=str2num(datal{1,1}{2,1});
            fa_left(sub,1)=str2num(datal{1,2}{33,1});

            qa_left(sub,1)=str2num(datal{1,2}{30,1});

            rd_left(sub,1)=str2num(datal{1,2}{36,1});

            ad_left(sub,1)=str2num(datal{1,2}{35,1});

            md_left(sub,1)=str2num(datal{1,2}{34,1});
        end


        for sub = 1:length(fnames_right)   
            fid = fopen(fnames_right{sub});
            datar = textscan(fid,'%s%s%s');
            fclose(fid);

            num_tracts_right(sub,1)=str2num(datar{1,1}{2,1});
            fa_right(sub,1)=str2num(datar{1,2}{33,1});

            qa_right(sub,1)=str2num(datar{1,2}{30,1});

            rd_right(sub,1)=str2num(datar{1,2}{36,1});

            ad_right(sub,1)=str2num(datar{1,2}{35,1});

            md_right(sub,1)=str2num(datar{1,2}{34,1});
        end
        
        % Let's get some descriptive stats going on here
        right_stats.(tracks{trk}) = [num_tracts_right,fa_right,qa_right,rd_right,ad_right,md_right];
        left_stats.(tracks{trk}) = [num_tracts_left,fa_left,qa_left,rd_left,ad_left,md_left];
        mean(right_stats.(tracks{trk}));
        mean(left_stats.(tracks{trk}));

        % Time to store data and plot it
        left_stats.(tracks{trk}) = array2table(left_stats.(tracks{trk}));
        left_stats.(tracks{trk}).Properties.VariableNames = {'num_tracts','fa','qa','rd','ad','md'};
        right_stats.(tracks{trk}) = array2table(right_stats.(tracks{trk}));
        right_stats.(tracks{trk}).Properties.VariableNames = {'num_tracts','fa','qa','rd','ad','md'};

        figure();
        tiledlayout(3,2);
        histogram(right_stats.(tracks{trk}).num_tracts); hold on; histogram(left_stats.(tracks{trk}).num_tracts);
        title('Number of tracts')
        nexttile
        histogram(right_stats.(tracks{trk}).fa); hold on; histogram(left_stats.(tracks{trk}).fa);
        title('Average FA')
        nexttile
        histogram(right_stats.(tracks{trk}).qa); hold on; histogram(left_stats.(tracks{trk}).qa);
        title('Average QA')
        nexttile
        histogram(right_stats.(tracks{trk}).rd); hold on; histogram(left_stats.(tracks{trk}).rd);
        title('Average RD')
        nexttile
        histogram(right_stats.(tracks{trk}).ad); hold on; histogram(left_stats.(tracks{trk}).ad);
        title('Average AD')
        nexttile
        histogram(right_stats.(tracks{trk}).md); hold on; histogram(left_stats.(tracks{trk}).md);
        title('Average MD')
    end
    if strcmp(tracks{trk},'Uncinate') == 0 
        fnames = filenames(fullfile(datadir,session,'results',tracks{trk},'*.stat.txt'));
        for sub = 1:length(fnames)
            if strcmp(fnames{trk}(1,1),'R') 
                fid = fopen(fnames{sub});
                data = textscan(fid,'%s%s%s');
                fclose(fid);

                num_tracts_left(sub,1)=str2num(datal{1,1}{2,1});
                fa_right(sub,1)=str2num(data{1,3}{15,1});

                qa_right(sub,1)=str2num(data{1,3}{9,1});

                rd_right(sub,1)=str2num(data{1,3}{21,1});

                ad_right(sub,1)=str2num(data{1,3}{19,1});

                md_right(sub,1)=str2num(data{1,3}{17,1});
                
                right_stats.(tracks{trk}) = [num_tracts_right,fa_right,qa_right,rd_right,ad_right,md_right];
                mean(right_stats.(tracks{trk}));
                
                right_stats.(tracks{trk}) = array2table(right_stats.(tracks{trk}));
                right_stats.(tracks{trk}).Properties.VariableNames = {'num_tracts','fa','qa','rd','ad','md'}; 
            else
                fid = fopen(fnames{sub});
                data = textscan(fid,'%s%s%s');
                fclose(fid);

                num_tracts_left(sub,1)=str2num(datal{1,1}{2,1});
                fa_left(sub,1)=str2num(data{1,3}{15,1});

                qa_left(sub,1)=str2num(data{1,3}{9,1});

                rd_left(sub,1)=str2num(data{1,3}{21,1});

                ad_left(sub,1)=str2num(data{1,3}{19,1});

                md_left(sub,1)=str2num(data{1,3}{17,1});
                
                left_stats.(tracks{trk}) = [num_tracts_left,fa_left,qa_left,rd_left,ad_left,md_left];
        
                mean(left_stats.(tracks{trk}));

                left_stats.(tracks{trk}) = array2table(left_stats.(tracks{trk}));
                left_stats.(tracks{trk}).Properties.VariableNames = {'num_tracts','fa','qa','rd','ad','md'};                
            end
            
        end
        
        if exist('right_stats') == 1
            curr_stats = right_stats;
        else
            curr_stats = left_stats;
        end
        
        figure();
        tiledlayout(3,2);
        title(strcat(tracks{trk},' track statistics'))
        histogram(curr_stats.(tracks{trk}).num_tracts);
        title('Number of tracts')
        nexttile
        histogram(curr_stats.(tracks{trk}).fa);
        title('Average FA')
        nexttile
        histogram(curr_stats.(tracks{trk}).qa);
        title('Average QA')
        nexttile
        histogram(curr_stats.(tracks{trk}).rd);
        title('Average RD')
        nexttile
        histogram(curr_stats.(tracks{trk}).ad);
        title('Average AD')
        nexttile
        histogram(curr_stats.(tracks{trk}).md);
        title('Average MD')  
    end

    clear fnames_left fnames_right left_stats right_stats num_tracts_right fa_right qa_right rd_right ad_right md_right num_tracts_left fa_left qa_left rd_left ad_left md_left
end
%%Script from Hartley lab members
%%edited by Gail Rosenbaum on 3/27/19
%%This script takes matlab output files from MBMF, computes relevant
%%behavioral variables for each subject (proportion of stays and switches
%%after common and rare transitions), creates individual-level and
%%group-level bar graphs of stay/switch probabilities, and creates an
%%output file that can be used to run a mixed model in R

clearvars

%%Directories containing output matlab files from MBMF task (should end in
%%_onsets.mat if the original code wasn't changed)
inputdir = '/Users/gailrosenbaum/Box Sync/MBMF/Scripts/BehavioralScripts/MATLAB/RawData/'; 

%where you want your output file to be saved, to submit to mixed model in R
outputdir ='/Users/gailrosenbaum/Box Sync/MBMF/Scripts/BehavioralScripts/MATLAB/RawData/Proc/';

files = dir(strcat(inputdir,'/*_onsets.mat'));

nGoodSubs = 0;
regDat=[];

for i=1:length(files)
    load([inputdir files(i).name]); %load subject data
    goodSub=1;
    
    %%%%%%%%%%%%% EXCLUSION CRITERIA %%%%%%%%%%%%%
    %Note that these criteria have not been consistently applied across
    %2-step studies. You may or may not want to implement them, but they
    %should help give you an idea of your data quality
    
    % does not currently exclude first n trials (although the stay-switch analysis below excludes 9 trials, as starttrial=10)
    
    % remove subjects with more than 15 missed trials   
    if sum(choice2==0) > 15 %%%Changed from 20, which was the initial exclusion criterion 
        fprintf('BAD_SUB: %d trials missed: %s \n', sum(choice2==0), name);
        goodSub=0;
    end
    % remove "uneven responders": prefer one stim over another in each
    % state,choosing unpreferred less than "tolerance" = .1
    tolerance = .1;
    if abs((sum(choice1==1)/sum(choice1~=0))-.5) > (.5 - tolerance)
        fprintf('BAD_SUB: Uneven S1 chooser: %s \n', name);
        goodSub=1; % if set to 1, turned off
    end
    if abs((length(intersect(find(state==2),find(choice2==1)))/length(intersect(find(state==2),find(choice2~=0))))-.5) > (.5 - tolerance) 
        fprintf('BAD_SUB: Uneven S2_1 chooser: %s \n', name);
        goodSub=1; % if set to 1, turned off
    end
    if abs((length(intersect(find(state==3),find(choice2==1)))/length(intersect(find(state==3),find(choice2~=0))))-.5) > (.5 - tolerance) 
        fprintf('BAD_SUB: Uneven S2_2 chooser: %s \n', name);
        goodSub=1; % if set to 1, turned off
    end

    % remove bad second stage choosers - if you choose a second stage stim
    % which rewarded you last time in the state less than 50% of the
    % time

    % below, if you stay more for common unrewarded than rewarded 

    
    if goodSub ==1
        common = double((choice1==1&state==2)|(choice1==2&state==3));
        common(common==0)=-1;
        rare= double((choice1==1&state==3)|(choice1==2&state==2));
        rare(rare==0)=-1;
        reward = money;
        reward(reward==0)=-1;
        validtrials = find(choice2~=0); %exclude missed trials from analysis
        common = common(validtrials); %include valid trials only
        rare = rare(validtrials);
        reward = reward(validtrials);
        stage1choice = choice1(validtrials);
        s1rts=rts1(validtrials);
        s2rts=rts2(validtrials);

        starttrial = 10; %trial at which to begin analysis
        r_c=[]; %did the participant stay after a rewarded trial that made a common transition? predict stay for all
        r_r=[]; %did the participant stay after a rare transition that was rewarded?; predict switch if MB, stay if MF
        u_c=[]; %predict switch for all
        u_r=[];  %predict switch if MB, stay if MF
        rcRT=[];
        rrRT=[];
        ucRT=[];
        urRT=[];
        % sort rts and stay choices into r_c, r_r, u_c, u_r
        for j = starttrial:length(validtrials)
            stay(j) = stage1choice(j)==stage1choice(j-1);
            lastwin(j) = reward(j-1);
            lasttransrare(j) = rare(j-1);
     
            if (lastwin(j) == 1) %rewarded
                if(lasttransrare(j) == -1) %common
                    rcRT(length(rcRT)+1)=s1rts(j);
                    if (stay(j)==1)
                        r_c(length(r_c)+1) = 1;
                    else
                        r_c(length(r_c)+1) = 0;
                    end;
                else %rare
                    rrRT(length(rrRT)+1)=s1rts(j);
                    if (stay(j)==1)
                        r_r(length(r_r)+1) = 1;
                    else
                        r_r(length(r_r)+1) = 0;
                    end;
                end
            else %unrewarded
                if(lasttransrare(j) == -1) %common
                    ucRT(length(ucRT)+1)=s1rts(j);
                    if (stay(j)==1)
                        u_c(length(u_c)+1) = 1;
                    else %switch
                        u_c(length(u_c)+1) = 0;
                    end
                else %rare
                    urRT(length(urRT)+1)=s1rts(j);
                    if (stay(j)==1)
                        u_r(length(u_r)+1) = 1;
                    else
                        u_r(length(u_r)+1) = 0;
                    end
                end
            end
        end
        if 1==1%(mean(u_c)-mean(r_c))<.05 %hacky way i implemented the exclusion of U_C > R_C
            nGoodSubs = nGoodSubs + 1; %if they pass, include in data record
            stay(1:starttrial-1)=[];
            lastwin(1:starttrial-1)=[];
            lasttransrare(1:starttrial-1)=[];

            s1rts=s1rts(starttrial:end);
            s2rts=s2rts(starttrial:end);
            currtransrare=rare(starttrial:end);

            %compile regression data to regress.dat file
            if length(name)==3
                id = repmat(str2num(name(end-2:end)),length(stay),1); 
            else
                id = repmat(str2num(name(end-3:end)),length(stay),1);
            end
            currRegDat=dataset({id,'subj'},{lastwin','lastwin'},{lasttransrare','lasttransR'},...
                {double(stay)','stay'},{(s1rts)','s1rt'},{(s2rts)','s2rt'},{(currtransrare)','currtransR'});
            regDat = [regDat;currRegDat]; %add to full set


            rc_means(i)=mean(r_c);
            uc_means(i)=mean(u_c);
            rr_means(i)=mean(r_r);
            ur_means(i)=mean(u_r);
            
            rcRTmeans(i)=mean(rcRT);
            ucRTmeans(i)=mean(ucRT);
            rrRTmeans(i)=mean(rrRT);
            urRTmeans(i)=mean(urRT);
            

            %plot subject data
            figure(1);
            subplot(7,10,i) %%%Change subplot dimensions so that you have enough for all of your data. here, the dimensions are 3x5, so if you have more than 15 files, you'll need to add (e.g., subplot(4,5,i) );
            bar([1 3],[rc_means(i),uc_means(i)],.45)
            hold on;
            bar([2 4],[rr_means(i),ur_means(i)],.45,'r')
            ylim([.05 1]);
            legend(name);
            hold off;
        else %if (goodsub == 0)
            fprintf('BAD_SUB: UC>RC: %s \n', name);
        end
        clear lastwin;
        clear lasttransrare;
        clear stay;
        clear id;
    end 
end

%%%plot group data; uncomment when running whole group
groupRC=mean(rc_means(:));
groupUC=mean(uc_means(:));
groupRR=mean(rr_means(:));
groupUR=mean(ur_means(:));
semRC=std(rc_means)/sqrt(length(rc_means));
semUC=std(uc_means)/sqrt(length(uc_means));
semRR=std(rr_means)/sqrt(length(rr_means));
semUR=std(ur_means)/sqrt(length(ur_means));

figure();
bar(1,groupRC,.8);
hold on;
errorbar(1,groupRC,semRC,'k');
bar(2,groupRR,.8,'r');
errorbar(2,groupRR,semRR,'k');
bar(3,groupUC,.8);
errorbar(3,groupUC,semUC,'k');
bar(4,groupUR,.8,'r');
errorbar(4,groupUR,semUR,'k');
ylim([.05 1]);
xlabel('reward''no reward')
hold off;

%%save MLM data file
export(regDat,'File',[outputdir '/AllSubsRegress.dat'],'delimiter',',');


nGoodSubs
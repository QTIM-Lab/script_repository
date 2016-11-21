function dcerecon_popAIF(input4Dfile, outputpath, maskfile)

% addpath /autofs/cluster/qtim/users/you2/DCEEstimation/code/
% addpath /autofs/cluster/mi2b2/4/you2/Utilities/NIfTI_matlab_20140122/
input4Dfile
if (nargin==2)
	maskfile='';
end

% load dependent parameters
DCE_OPTS = [  0    1    0    1   0   1   1   2    1    1   1    1    1    0 ];
T1_for_fit = 1000; %assumed T1 in tissue
T1_fixed = 1000; %in milliseconds - NOTE THAT THIS MAY BE USED FOR THE CONCENTRATION CALCULATION INSTEAD OF T1 MAPS
r1_Cagent = 0.0039; %Relaxometry of Contrast Agent used jkc

TR  = 6.8; %5; % Repetition time in millisec!
alpha=30;alpha_rad = (pi/180).*alpha;
Total_scan_time_mins = 6; %11; %jkc
%nr_of_frames =  Total_scan_time_mins/FR_mins+1;
firstbaseline = 1;
firstBaseline = firstbaseline;



                if (exist(input4Dfile))
                        dce4D=load_untouch_nii(input4Dfile);
% 			mask=load_untouch_nii(maskfile);
                        xsize=size(dce4D.img,1);
                        ysize=size(dce4D.img,2);
                        zsize=size(dce4D.img,3);
                        tsize=size(dce4D.img,4);
			tsize
			%if (tsize==200)
			%	lastbaseline = 24; % for TIV_01/VISIT_01 and VISIT_02
			%	FR=2.5;
			%end
% 			if (tsize==250)
% 				lastbaseline = 42; % all others
% 				FR=1.6;
% 			end
            if (tsize==60)
                lastbaseline = 28;
                FR=6;
            else
                'wrong tsize'
                return
            end
            
			lastBaseline = lastbaseline;
			lpbs=lastbaseline+1;     %dynamic #
			FR_mins = FR/60;         % in minutes


                        ktransmap=zeros(xsize,ysize,zsize,tsize);
                        vemap=zeros(xsize,ysize,zsize,tsize);
                        aucmap=zeros(xsize,ysize,zsize,tsize);


                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        % generate aif (popAIF)
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        [AIF]=generateAIF(tsize,FR,lpbs);


                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        % calculate ktrans, ve at each voxel
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        for z=1:zsize
                          fprintf('slice %d\n', z);
                          for x=1:xsize
                            for y=1:ysize
%                                                         for z=1:zsize
%                           fprintf('slice %d\n', z);
%                           for x=40:42
%                             for y=50:52
				if (1 == 1)
%                 if (mask.img(x,y,z)>0)
                                % signal
                                tic
                                s=dce4D.img(x,y,z,:);
				input_signal_4D=zeros(1,1,tsize);
                                for n=1:tsize
                                        input_signal_4D(1,1,n)=s(1,1,1,n);
                                end
%                                 input_signal_4D
                                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                % convert 4D signal to 4D concentration (reference to
                                % Step2a_DRO_signal_to_concentration.m)
                                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                signal_4D = input_signal_4D;
                                baselineVol_3D = mean(signal_4D(:,:,firstBaseline:lastBaseline),3);
                                if (baselineVol_3D == 0)
                                ktransmap(x,y,z)=0;
                                vemap(x,y,z)=0;
                                aucmap(x,y,z)=0;
                                continue
                                end
                                relSignal_4D = zeros(size(signal_4D));
                                for i=1:size(signal_4D,3)
                                    relSignal_4D(:,:,i) = signal_4D(:,:,i)./baselineVol_3D;
                                end
%                                 relSignal_4D
                                R1pre = 1./T1_for_fit;   %1/msec
                                a = exp(-TR.*R1pre);
                                TERM = (1-a)./(1-a.*cos(alpha_rad)); 
%                                 TERM
                                y_4D = relSignal_4D.*(repmat(TERM,[size(signal_4D,1),size(signal_4D,2),size(signal_4D,3)]));
                                % Use y_4D to calculate CA concentration:
%                                 y_4D
                                gd_conc_4D = zeros(size(relSignal_4D));
                                for i=1:size(relSignal_4D,3);
                                    y_3D = squeeze(y_4D(:,:,i));
                                    gd_log_term = (y_3D-1)./(a.*(y_3D.*cos(alpha_rad)-1));
                                    gd_conc_4D(:,:,i) = -(1./(r1_Cagent*TR)) .* log(gd_log_term);
                                end
%                                 gd_conc_4D
                                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                % fit ktrans and Ve (simplex)
                                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                obs_conc=squeeze(gd_conc_4D);
%                                 init_params=[1 1];
                                init_params=[-2, 0.1];
                                [kin_par est_conc] = Step4b_Simplex_Fit_Kinetic_Tofts(obs_conc, AIF, init_params, FR_mins);
                                k1=exp(kin_par(1)); ktrans=k1;
                                Ve=1/(1+exp(-kin_par(2)));
                                k2 = k1/Ve;
                                toc
                                fprintf('at (%d, %d, %d), Ve=%f, ktrans=%f\n', x, y, z, Ve, ktrans);
                                ktransmap(x,y,z)=ktrans;
                                vemap(x,y,z)=Ve;
                                aucmap(x,y,z)=trapz(obs_conc)/trapz(AIF);
				end
                            end
                          end
                        end

                        % save ktrans, ve, auc maps
                        tmp=dce4D;
                        tmp.hdr.dime.dim(1)=3;
                        tmp.hdr.dime.dim(5)=1;
			tmp.hdr.dime.pixdim(1)=1;
			tmp.hdr.dime.datatype=16;
			tmp.hdr.dime.bitpix=32;  % make sure it is a float image
			tmp.hdr.dime.cal_max=0;
			tmp.hdr.dime.glmax=0;
            
                        path = input4Dfile
                        visit = path(end-31:end-24)
                        patient = path(end-38:end-33)
                     
                        
                        tmp.img = tmp.img(:,:,:,1)
                        fn=['dce_01.nii.gz']
                        save_untouch_nii(tmp,strcat(patient, '_', visit, '_', fn));

                        tmp.img=ktransmap;
                        fn=['ktrans.nii.gz']
                        save_untouch_nii(tmp,strcat(patient, '_', visit, '_', fn));

                        tmp.img=vemap;
                        fn=['ve.nii.gz']
                        save_untouch_nii(tmp,strcat(patient, '_', visit, '_', fn));

                        tmp.img=aucmap;
                        fn=['auc.nii.gz']
                        save_untouch_nii(tmp,strcat(patient, '_', visit, '_', fn));
                                     
		end                     

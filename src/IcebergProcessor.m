classdef IcebergProcessor
    properties
        signal
        level
        angle
        IR
        configurationSetup
        VBAP_DSER_Part
        Amb_ERLR_Part
        
        iFs
        signalAmb
        iChannel
    end
    
    methods
        function obj = IcebergProcessor(signal, level, angle, IR, configurationSetup)
            if nargin < 5 || isempty(configurationSetup)
                configurationSetup = obj.getDefaultConfigurationSetup();
            end
            obj.signal = signal;
            obj.level = level;
            obj.angle = angle;
            obj.IR = IR;
            obj.configurationSetup = configurationSetup;
            
            % Initialize itaAudio objects
            obj.VBAP_DSER_Part = itaAudio();
            obj.Amb_ERLR_Part = itaAudio();
            obj.VBAP_DSER_Part.samplingRate = signal.samplingRate;
            obj.Amb_ERLR_Part.samplingRate = signal.samplingRate;
        end
        
        function final_audio = process(obj)
            % Array Settings
            obj.configurationSetup.activeLSNumbers = zeros(1, length(obj.configurationSetup.ls_dir));
            for indexx = 1:length(obj.configurationSetup.ls_dir)
                iArrayP = find(obj.configurationSetup.LSArray == obj.configurationSetup.ls_dir(indexx, 1), 1);
                obj.configurationSetup.activeLSNumbers(indexx) = iArrayP;
            end
            
            obj.iFs = obj.signal.samplingRate; % Sample Frequency
            obj.signal = ita_normalize_dat(obj.signal, 'allchannels', 'true');
            
            %% Ambisonics part
            %Load the most recent calibration
            obj.configurationSetup.Level_Factor = obj.configurationSetup.new_Level_Factor;
            obj.configurationSetup.lsdBperVolt = (20*log10((obj.configurationSetup.iFactor)/2e-5));
            if obj.level > 95
                error ('Level max is 95 dB');
            end
            signal_to_play = obj.signal;
            
            obj.configurationSetup.activeLSNumbers = zeros(1,length(obj.configurationSetup.ls_dir));
            for indexActiveLSNumbers= 1:length(obj.configurationSetup.ls_dir)
                iArrayP = find(obj.configurationSetup.LSArray==obj.configurationSetup.ls_dir(indexActiveLSNumbers,1),1);
                obj.configurationSetup.activeLSNumbers(indexActiveLSNumbers) = iArrayP;
            end
            for iCount = obj.configurationSetup.activeLSNumbers
                Interpolation(:,iCount) = pchip(obj.configurationSetup.iLoudspeakerFreqFilter(iCount).freqVector,...
                    obj.configurationSetup.iLoudspeakerFreqFilter(iCount).freq,signal_to_play.freqVector);
            end
            
            frequencyFilter = itaAudio(Interpolation,obj.iFs,'freq');
            new_page = zeros(length(signal_to_play.time),1);
            stimulus = zeros(length(signal_to_play.time),1);
            signal_run = zeros(length(signal_to_play.time),1);
            
            %% Select the filter fo the audio according to the LS (Virtual Loudspeakers are filtered with the Nearest Speaker NSP)
            % This is to fit level the ambisonics audio part
            if obj.angle > 90 && obj.angle <= 180
                if  obj.angle <= 135
                    obj.iChannel = obj.configurationSetup.activeLSNumbers(find(obj.configurationSetup.ls_dir==90,1));
                else
                    obj.iChannel = obj.configurationSetup.activeLSNumbers(find(obj.configurationSetup.ls_dir==180,1));
                end
            elseif obj.angle > 180 && obj.angle <= 270
                if  obj.angle <= 225
                    obj.iChannel = obj.configurationSetup.activeLSNumbers(find(obj.configurationSetup.ls_dir==180,1));
                else
                    obj.iChannel = obj.configurationSetup.activeLSNumbers(find(obj.configurationSetup.ls_dir==270,1));
                end
            elseif obj.angle > 270 && obj.angle <= 360
                if  obj.angle <= 315
                    obj.iChannel = obj.configurationSetup.activeLSNumbers(find(obj.configurationSetup.ls_dir==270,1));
                else
                    obj.iChannel = obj.configurationSetup.activeLSNumbers(find(obj.configurationSetup.ls_dir==0,1));
                end
            elseif obj.angle >= 0 && obj.angle <= 90
                if  obj.angle <= 45
                    obj.iChannel = obj.configurationSetup.activeLSNumbers(find(obj.configurationSetup.ls_dir==0,1));
                else
                    obj.iChannel = obj.configurationSetup.activeLSNumbers(find(obj.configurationSetup.ls_dir==90,1));
                end
            end
            
            % Calculate the absolute difference between the user's number and each element of the vector
            diff_vec = abs(obj.configurationSetup.ls_dir(1,:) - obj.angle);
            
            % Find the index of the closest element
            [~, obj.iChannel] = min(diff_vec);
            
            signal_run_ita_FILTER       = itaAudio(zeros(length(signal_to_play.time),1),obj.iFs,'time');
            signal_run_ita_FILTER_LEVEL = itaAudio(zeros(length(signal_to_play.time),1),obj.iFs,'time');
            
            if isnan(signal_to_play.time)
                new_page(:,1) = 0;
            else
                new_page(:,1) = signal_to_play.time;
            end
            scaler = (sqrt(mean(new_page(:,1).^2)))*2;
            stimulus(:,1) = new_page(:,1)./scaler;
            if obj.level ~= 'n' %Verify this
                signal_run(:,1) = stimulus(:,1).*repmat(10.^((obj.level - obj.configurationSetup.lsdBperVolt )./20),length(new_page(:,1)),1);
            else
                signal_run(:,1)= signal_to_play.time;
            end
            signal_run_ita_dB = itaAudio(signal_run,obj.iFs,'time');
            filtered = ita_multiply_spk(signal_run_ita_dB,frequencyFilter.ch(obj.iChannel));
            signal_run_ita_FILTER.time(:,1) = filtered.time;
            signal_run_ita_FILTER_LEVEL.time(:,1) = signal_run_ita_FILTER.time.*(obj.configurationSetup.Level_Factor(obj.iChannel));
            
            
            obj.signalAmb = signal_run_ita_FILTER_LEVEL;
            
            %% end Ambisonics normalization part (Check all this - smells fishy)
            
            % Rest of the code...
            [D4,~] = ambiDecoder(obj.configurationSetup.ls_dir,'SAD',1,1);                                 % Create Decoder Matrix with n=4 LS SAD decoder!
            omnichannelIR = ita_split(obj.IR,1);                                %Select/get Omnidirectional IR
            %% Center time
            [IR_Early, shiftIndex] = ita_time_shift(omnichannelIR,'auto'); % Pull
            %
            %             if iIR == 3
            %                 IR_Early = ita_time_window(IR_Early,[0 .01],'time','windowType','rectwin'); %It does not make sense the cTime Only DS from Odeon simulation
            %                 DSER = ita_time_shift(IR_Early,abs(shiftIndex));            % Push
            %                 cTime = 0.01;
            %             else
            centerTime = ita_roomacoustics(omnichannelIR,'Center_Time','broadbandAnalysis',1);
            cTime = centerTime.Center_Time.freq;
            if isnan(cTime)
                centerTime = ita_roomacoustics(ita_normalize_dat(omnichannelIR),'Center_Time','broadbandAnalysis',1,'edcMethod','noCut');
                cTime = centerTime.Center_Time.freq;
            end
            IR_Early = ita_time_window(IR_Early,[0 cTime],'time','windowType','hann');
            % Shift back to fit the future composition
            DSER = ita_time_shift(IR_Early,abs(shiftIndex));            % Push
            %             end
            
            final_audio = ita_add(obj.VBAP_DSER_Part, obj.Amb_ERLR_Part);
            
            if length(obj.configurationSetup.LSArray) > final_audio.dimensions
                final_audio.time(:, final_audio.dimensions+1:length(obj.configurationSetup.LSArray)) = zeros(final_audio.nSamples, length(obj.configurationSetup.LSArray)-(final_audio.dimensions));
            end
            
            %% Signal with Omnidirectional reverberation
            convolved_DSER_signal = ita_convolve(obj.signal,DSER);
            %% Calibrated VBAP to specified Sound Pressure Level
            %             VBAP_DS =   set_level_vbap_fly_in(convolved_DSER_signal,level,angle,configurationSetup);
            
            
            %% Load values from the most recent calibration
            Level_Factor = obj.configurationSetup.new_Level_Factor; %LS Level adjust
            obj.iFs = obj.signal.samplingRate;                          %Sample Frequency
            lsdBperVolt = (20*log10((obj.configurationSetup.iFactor)/2e-5)); %Factor
            %% Prepare Signal
            signal_to_play = itaAudio();        %Empty itaAudio object
            signal_to_play.samplingRate = obj.iFs;  %set FS (if you use other than 44.1 this step is essential)
            signal_to_play.fftDegree = log2(obj.iFs*obj.signal.trackLength); %At that point you got the FFTdegree to seconds already, right? t = (2^(fftDegree))/fs
            for nSignalsIndex = 1:obj.signal.nChannels
                signal_vector(nSignalsIndex,:,:) = obj.ring_VBAP(obj.iFs,obj.signal.time(:,nSignalsIndex),obj.angle(nSignalsIndex),obj.configurationSetup);
                signal_ita(nSignalsIndex)  = itaAudio(squeeze(signal_vector(nSignalsIndex,:,:)),obj.iFs,'time');
                signal_to_play = ita_add(signal_to_play,signal_ita(nSignalsIndex));
            end
            
            %%
            activeLSNumbers = zeros(1,length(obj.configurationSetup.ls_dir));
            for indexx= 1:length(obj.configurationSetup.ls_dir)
                iArrayP = find(obj.configurationSetup.LSArray==obj.configurationSetup.ls_dir(indexx,1),1);
                activeLSNumbers(indexx) = iArrayP;
            end
            for iCount = activeLSNumbers
                Interpolation(:,iCount) = pchip(obj.configurationSetup.iLoudspeakerFreqFilter(iCount).freqVector,...
                    obj.configurationSetup.iLoudspeakerFreqFilter(iCount).freq,signal_to_play.freqVector);
            end
            
            frequencyFilter = itaAudio(Interpolation,obj.iFs,'freq');
            
            %%
            new_page = zeros(length(signal_to_play.time),length(activeLSNumbers));
            scaler = zeros(1,length(activeLSNumbers));
            stimulus = zeros(length(signal_to_play.time),length(activeLSNumbers));
            signal_run = zeros(length(signal_to_play.time),length(activeLSNumbers));
            levels = zeros(1,max(activeLSNumbers));
            
            %% VBAP balance receive the correct SPL from the two speakers
            
            angle_space = 360/length(obj.configurationSetup.ls_dir); %
            %Clockwise adjust to match odeon
            
            if obj.angle > 90 && obj.angle<= 180 %Desired presentation to sequential LS
                if  obj.angle<= 135
                    %     s1 = 1; s2 = 19;
                    s1 = activeLSNumbers(obj.configurationSetup.ls_dir(:,1)==180);
                    s2 = activeLSNumbers(obj.configurationSetup.ls_dir(:,1)==90);
                else
                    s1 = activeLSNumbers(obj.configurationSetup.ls_dir(:,1)==90);
                    s2 = activeLSNumbers(obj.configurationSetup.ls_dir(:,1)==180);
                end
            elseif obj.angle> 180 && obj.angle<= 270
                if  obj.angle<= 225
                    s1 = activeLSNumbers(obj.configurationSetup.ls_dir(:,1)==180);
                    s2 = activeLSNumbers(obj.configurationSetup.ls_dir(:,1)==270);
                else
                    s1 = activeLSNumbers(obj.configurationSetup.ls_dir(:,1)==180);
                    s2 = activeLSNumbers(obj.configurationSetup.ls_dir(:,1)==270);
                end
            elseif obj.angle> 270 && obj.angle<= 360
                if  obj.angle<= 315
                    s1 = activeLSNumbers(obj.configurationSetup.ls_dir(:,1)==270);
                    s2 = activeLSNumbers(obj.configurationSetup.ls_dir(:,1)==0);
                else
                    s1 = activeLSNumbers(obj.configurationSetup.ls_dir(:,1)==0);
                    s2 = activeLSNumbers(obj.configurationSetup.ls_dir(:,1)==270);
                end
            elseif (obj.angle>= 0 && obj.angle<= 90)
                if  obj.angle<= 45
                    s1 = activeLSNumbers(obj.configurationSetup.ls_dir(:,1)==0);
                    s2 = activeLSNumbers(obj.configurationSetup.ls_dir(:,1)==90);
                else
                    s1 = activeLSNumbers(obj.configurationSetup.ls_dir(:,1)==90);
                    s2 = activeLSNumbers(obj.configurationSetup.ls_dir(:,1)==0);
                end
            end
            
            %determine panning levels
            %Coherent sum (6 dB) , so ^2 to calculate the correct energy level factor %Otherwise just sin/cos to pan law
            s1_level = (sin(abs(rem(obj.angle,angle_space)/angle_space*(pi/2))))^2;
            s2_level = (cos(abs(rem(obj.angle,angle_space)/angle_space*(pi/2))))^2;
            
            [s1_level, s2_level] = deal(max([s1_level s2_level]),min([s1_level  s2_level]));
            if obj.level ~= 'n'
                for idx = 1:length(obj.level)
                    %         [s1(idx), s2(idx) s1_level(idx) s2_level(idx)] = equal_power_pan(obj.angle(idx));
                    if isinf(20*log10(10^((obj.level(idx)/20))*s1_level(idx)))==1
                        levels(s1(idx)) = 0;
                    else
                        levels(s1(idx)) = 20*log10(10^((obj.level(idx)/20))*s1_level(idx));
                    end
                    if isinf(20*log10(10^((obj.level(idx)/20))*s2_level(idx)))
                        levels(s2(idx)) = 0;
                    else
                        levels(s2(idx)) = 20*log10(10^((obj.level(idx)/20))*s2_level(idx));
                    end
                end
            end
            signal_run_ita_FILTER       = itaAudio(zeros(length(signal_to_play.time),length(activeLSNumbers)),obj.iFs,'time');
            signal_run_ita_FILTER_LEVEL = itaAudio(zeros(length(signal_to_play.time),length(activeLSNumbers)),obj.iFs,'time');
            
            %% What is happening here?
            for idx = 1:length(activeLSNumbers)
                if isnan(signal_to_play.ch(idx).time)==1
                    new_page(:,idx) = 0;
                    fprintf('NaN Channel, take a look ch: %i\n',idx);
                else
                    
                    new_page(:,idx) = signal_to_play.ch(idx).time;
                end
                scaler(idx) = (sqrt(mean(new_page(:,idx).^2)))*2;
                if (scaler(idx)~=0) %Avoid NaN
                    stimulus(:,idx) = new_page(:,idx)./scaler(idx);
                end
                %     if level ~= 'n'
                signal_run(:,idx) = stimulus(:,idx).*repmat(10.^((levels(activeLSNumbers(idx)) - lsdBperVolt )./20),length(new_page(:,idx)),1);
                %     else
                %         signal_run(:,idx)= signal_to_play.ch(idx).time;
                %     end
                signal_run_ita_dB = itaAudio(signal_run,obj.iFs,'time');
                filtered = ita_multiply_spk(signal_run_ita_dB.ch(idx),frequencyFilter.ch(activeLSNumbers(idx)));
                signal_run_ita_FILTER.time(:,idx) = filtered.time;
                signal_run_ita_FILTER_LEVEL.time(:,idx) = signal_run_ita_FILTER.ch(idx).time.*(Level_Factor(activeLSNumbers(idx)));
            end
            %%
            VBAP_DS = signal_run_ita_FILTER_LEVEL;
            
            
            VBAP_DS =   VBAP_DS*max(DSER.time);
            %% fetch VBAP Direct Sound to the Array
            
            for i = 1:length(activeLSNumbers)
                obj.VBAP_DSER_Part.time(:,activeLSNumbers(i)) = VBAP_DS.time(:,i);
            end
            %% Shift IR to the begining
            [IR_Late, shiftIndex] = ita_time_shift(obj.IR,'auto');
            % Crop out the the Direc Sound+Early Reflections (Center Time) (May work
            % with window as well, but the energy balance need to be verified in this
            % case
            IR_Late = ita_time_crop(IR_Late,[cTime 0],'time');
            % Shift back to fit the future composition
            resyncSamples = obj.IR.nSamples-IR_Late.nSamples;               % Timefactor
            IR_Late.time = [zeros(resyncSamples,4); IR_Late.time];      % Adjusting time
            IR_Late = ita_time_window(IR_Late,...                       %Get rid of non linearities
                [0 (IR_Late.trackLength-0.05)],'time',...  %Get rid of non linearities
                'windowType','rectwin');                   %Get rid of non linearities
            IR_Late = ita_time_shift(IR_Late,abs(shiftIndex));          % Push
            
            %% Ambisonics
            sinal_Ambisonics_ERLR   = ita_convolve(obj.signalAmb,IR_Late);         % Signal with Late Reverberation
            %% Calibrated Ambisonics to specified Sound Pressure Level
            %Decode to Loudspeaker array
            Ambisonics_ERLRSignal = itaAudio(decodeBformat(...
                sinal_Ambisonics_ERLR.time,...
                D4),obj.iFs,'time');
            % fetch Ambisonics Late reverberation to the Array
            for i = 1:length(activeLSNumbers)
                obj.Amb_ERLR_Part.time(:,activeLSNumbers(i))  = Ambisonics_ERLRSignal.time(:,i);
            end
            %% Combine DS ER and LR
            final_audio = ita_add(obj.VBAP_DSER_Part,obj.Amb_ERLR_Part);    %%
            if length(obj.configurationSetup.LSArray) >final_audio.dimensions
                final_audio.time(:,final_audio.dimensions+1:length(obj.configurationSetup.LSArray)) = zeros(final_audio.nSamples,length(obj.configurationSetup.LSArray)-(final_audio.dimensions));
            end
            
        end

        function [pansig,padsig,sig] = ring_VBAP(obj,fs,iSignal,iAngle,configurationSetup)
            %% VBAP Archontis
            % initial parameters (Here we can even move the sound source... I will
            % update this to have as an option, now is static) .Sergio.
            blocksize = fs/18.3750; % (~50msec)
            if fs == 48000
                blocksize = fs/20; % (~50msec)
            elseif fs == 96000
                blocksize = fs/40; % (~50msec)
            end
            hopsize = blocksize/2; % panning hopsize (update the panning values twice per bocksize)
            ls_num = length(configurationSetup.ls_dir);
            sig = iSignal;
            %% Parametros do sinal
            Lsig = length(sig);
            Nhop = ceil(Lsig/hopsize) + 2;
            padsig = [zeros(hopsize,1); sig ;zeros(Nhop*hopsize - Lsig - hopsize,1)]; % zero padding
            pansig = zeros(size(padsig,1), ls_num);
            %% define the trajectory for panning (Quando sinal esta girando girando...)
            static = ones(length((0:(Nhop-1)-1)'*(9*360)/(Nhop-1)),1); %verify here
            iAngleCount = iAngle;
            azis = iAngleCount*static;
            eles = 0*azis;
            % precompute VBAP matrix inversion (I LOVE this part)
            ls_groups = findLsPairs(configurationSetup.ls_dir(:,1)); %it can return also a mesh for plotting
            layoutInvMtx = invertLsMtx(configurationSetup.ls_dir(:,1), ls_groups);
            %% do panning of signal to defined trajectory with overlap-add method
            counter = 1;
            window = hanning(blocksize);
            spread = 0;
            
            for idx = 0:hopsize:(Nhop-2)*hopsize
                winsig = padsig(idx+(1:blocksize),1).*window;
                azi = azis(counter);
                ele = eles(counter);
                gains   = vbap([azi ele], ls_groups, layoutInvMtx, spread); %Using VBAP
                gains2  = vbip([azi ele], ls_groups, layoutInvMtx, spread); %Using VBIP
                panwinsig = winsig*gains;
                pansig(idx+(1:blocksize),:) = pansig(idx+(1:blocksize),:) + panwinsig;
                counter = counter+1;
            end
            % truncate loudspeaker signals to original length (omit zeropadding)
            pansig = pansig(hopsize+(1:Lsig),:);
            
            
        end
    end
    
    

    methods (Access = private)


        
        function configurationSetup = getDefaultConfigurationSetup(obj)
            % Define and return the default configurationSetup
            
            
            configurationSetup.Channel_1 =  180; %The back loudspeaker is connected to output number 1 of the sound card.
            configurationSetup.Channel_7 =  270; %The left loudspeaker is connected to output number 7 of the sound card.
            configurationSetup.Channel_13 = 000; %The front loudspeaker is connected to output number 13 of the sound card.
            configurationSetup.Channel_19 = 090; %The right loudspeaker is connected to output number 19 of the sound card.
            
            configurationSetup.TotalNumberOfLoudspeakers = 24;
            configurationSetup.new_Level_Factor = ones(1,24);
            configurationSetup.iFactor = 1;
            
            b = itaResult;
            b.freqData = ones(27,1);
            b.freqVector = [50;62.5;80;100;125;155;200;250;315;400;500;630;800;1000;1250;1600;2000;2500;3150;4000;5000;6350;8000;10000;12500;16000;20000];
            
            for iLSpeaker = 1:24
                configurationSetup.iLoudspeakerFreqFilter(iLSpeaker) = b;
            end
            clear b
            configurationSetup.spaceBetweenLS = 360/configurationSetup.TotalNumberOfLoudspeakers;
            
            
            % To be used in 4 LS Hybrid mode we need to specify the LS index numbers
            configurationSetup.ls_dir = [configurationSetup.Channel_1 configurationSetup.Channel_7 configurationSetup.Channel_13 configurationSetup.Channel_19; 0 0 0 0]';
            %Right LS is 90 Degrees / Left is 270
            configurationSetup.LSArray = [mod(configurationSetup.Channel_1:configurationSetup.spaceBetweenLS:360,360) configurationSetup.spaceBetweenLS:configurationSetup.spaceBetweenLS:(configurationSetup.Channel_1-configurationSetup.spaceBetweenLS)]'; %Clockwise configuration
            %Right LS is 270 Degrees / Left is 90 counter = clockwise configuration
            % configurationSetup.LSArray = [mod(Channel_1:-spaceBetweenLS:0,360) (Channel_1-spaceBetweenLS):-spaceBetweenLS:spaceBetweenLS]';
        end
    end
end

        

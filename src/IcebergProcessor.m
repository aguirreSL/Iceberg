classdef IcebergProcessor
    properties
        configurationSetup
    end

        properties (Access = private)
        VBAP_DSER
        Ambisonics_LR
        processed_audio
        setup_audio
    end
    
    methods
        function obj = IcebergProcessor(signal, level, angle, IR, configurationSetup)
            
            if level > 95
                error ('Level max is 95 dB');
            end
            
            if nargin < 5 || isempty(configurationSetup)
                configurationSetup = obj.getDefaultConfigurationSetup();
            end

            configurationSetup.signal = signal;
            configurationSetup.level = level;
            configurationSetup.IR = IR;
            
            if angle <= 360
            angle = round(abs(angle)/5) * 5;
            else
                disp('Max angle is 360')
            end
            configurationSetup.angle = angle;

            obj.configurationSetup = configurationSetup;
            % Initialize itaAudio objects
            obj.VBAP_DSER = itaAudio();
            obj.Ambisonics_LR = itaAudio();
            obj.VBAP_DSER.samplingRate = signal.samplingRate;
            obj.Ambisonics_LR.samplingRate = signal.samplingRate;
        end
        
        function [final_audio, setup_audio] = process(obj)
    
            signal = ita_normalize_dat(obj.configurationSetup.signal, 'allchannels', 'true');

            obj.configurationSetup.Level_Factor = obj.configurationSetup.new_Level_Factor;
            obj.configurationSetup.lsdBperVolt = (20*log10((obj.configurationSetup.iFactor)/2e-5));
            obj.configurationSetup.activeLSNumbers = zeros(1,length(obj.configurationSetup.ls_dir));

            for indexActiveLSNumbers= 1:length(obj.configurationSetup.ls_dir)
                iArrayP = find(obj.configurationSetup.LSArray==obj.configurationSetup.ls_dir(indexActiveLSNumbers,1),1);
                obj.configurationSetup.activeLSNumbers(indexActiveLSNumbers) = iArrayP;
            end
            if ~ obj.configurationSetup.default 
            for iCount = obj.configurationSetup.activeLSNumbers
                Interpolation(:,iCount) = pchip(obj.configurationSetup.iLoudspeakerFreqFilter(iCount).freqVector,...
                    obj.configurationSetup.iLoudspeakerFreqFilter(iCount).freq,signal.freqVector);
            end
            frequencyFilter = itaAudio(Interpolation,signal.samplingRate,'freq');
            end
            diff_vec = abs(obj.configurationSetup.ls_dir(1,:) - obj.configurationSetup.angle);
            
            % Find the index of the closest element
            [~, iChannel] = min(diff_vec);
            
            new_page = zeros(length(signal.time),1);
            stimulus = zeros(length(signal.time),1);
            signal_run = zeros(length(signal.time),1);
            
            signal_run_ita_FILTER       = itaAudio(zeros(length(signal.time),1),signal.samplingRate,'time');
            signal_run_ita_FILTER_LEVEL = itaAudio(zeros(length(signal.time),1),signal.samplingRate,'time');
            
            if isnan(signal.time)
                new_page(:,1) = 0;
            else
                new_page(:,1) = signal.time;
            end

            scaler = (sqrt(mean(new_page(:,1).^2)))*2;
            stimulus(:,1) = new_page(:,1)./scaler;

            if obj.configurationSetup.level ~= 'n' %Verify this
                signal_run(:,1) = stimulus(:,1).*repmat(10.^((obj.configurationSetup.level - obj.configurationSetup.lsdBperVolt )./20),length(new_page(:,1)),1);
            else
                signal_run(:,1)= signal.time;
            end

            signal_run_ita_dB = itaAudio(signal_run,signal.samplingRate,'time');
            
            if obj.configurationSetup.default 
                signal_run_ita_FILTER.time(:,1) = signal_run_ita_dB.time;
            else
            filtered = ita_multiply_spk(signal_run_ita_dB,frequencyFilter.ch(iChannel));
            signal_run_ita_FILTER.time(:,1) = filtered.time;
            end
            signal_run_ita_FILTER_LEVEL.time(:,1) = signal_run_ita_FILTER.time.*(obj.configurationSetup.Level_Factor(iChannel));
            scaledSignal = signal_run_ita_FILTER_LEVEL;
                        
            [D4,~] = ambiDecoder(obj.configurationSetup.ls_dir,'SAD',1,1);                                 % Create Decoder Matrix with n=4 LS SAD decoder!
            
            omnichannelIR = ita_split(obj.configurationSetup.IR,1); 
            centerTime = ita_roomacoustics(omnichannelIR,'Center_Time','broadbandAnalysis',1);
            cTime = centerTime.Center_Time.freq;
            %if iIR zero RT
            %IR_Early = ita_time_window(IR_Early,[0 .01],'time','windowType','rectwin'); %It does not make sense the cTime Only DS from Odeon simulation
            %DSER = ita_time_shift(IR_Early,abs(shiftIndex));            % Push
            %cTime = 0.01; else--->>>
            if isnan(cTime)
                centerTime = ita_roomacoustics(ita_normalize_dat(omnichannelIR),'Center_Time','broadbandAnalysis',1,'edcMethod','noCut');
                cTime = centerTime.Center_Time.freq;
            end
            %% Shift IR to the begining
            [IR_Late, shiftIndex] = ita_time_shift(obj.configurationSetup.IR,'auto');
            % Crop out the the Direc Sound+Early Reflections (Center Time) (May work
            % with window as well, but the energy balance need to be verified in this
            % case
            IR_Late = ita_time_crop(IR_Late,[cTime 0],'time');
            
            % Shift back to fit the future composition
            resyncSamples = obj.configurationSetup.IR.nSamples-IR_Late.nSamples;           % Timefactor
            IR_Late.time = [zeros(resyncSamples,4); IR_Late.time];      % Adjusting time
            IR_Late = ita_time_window(IR_Late, [0 (IR_Late.trackLength-0.05)],'time','windowType','rectwin'); %Get rid of non linearities
            IR_Late = ita_time_shift(IR_Late,abs(shiftIndex));          % Push
            
            %% Ambisonics
            sinal_Ambisonics_ERLR   = ita_convolve(scaledSignal,IR_Late);         % Signal with Late Reverberation
            Ambisonics_LRSignal   = itaAudio(decodeBformat(sinal_Ambisonics_ERLR.time,D4),signal.samplingRate,'time');
            
            for i = 1:length(obj.configurationSetup.activeLSNumbers)
                obj.Ambisonics_LR.time(:,obj.configurationSetup.activeLSNumbers(i))  = Ambisonics_LRSignal.time(:,i);
            end

            [IR_Early, shiftIndex] = ita_time_shift(omnichannelIR,'auto'); % Pull
            IR_Early = ita_time_window(IR_Early,[0 cTime],'time','windowType','hann');
            DSER = ita_time_shift(IR_Early,abs(shiftIndex));            % Push
            convolved_DSER_signal = ita_convolve(signal,DSER);

            %% Prepare Signal
            signal_to_play = itaAudio();                     
            signal_to_play.samplingRate = signal.samplingRate;           
            signal_to_play.fftDegree = log2(signal.samplingRate*signal.trackLength);
            for nSignalsIndex = 1:signal.nChannels
                signal_vector(nSignalsIndex,:,:) = obj.ring_VBAP(signal.samplingRate,signal.time(:,nSignalsIndex),angle(nSignalsIndex),obj.configurationSetup);
                signal_ita(nSignalsIndex)  = itaAudio(squeeze(signal_vector(nSignalsIndex,:,:)),signal.samplingRate,'time');
                signal_to_play = ita_add(signal_to_play,signal_ita(nSignalsIndex));
            end
            
            new_page    = zeros(length(signal_to_play.time),length(obj.configurationSetup.activeLSNumbers));
            scaler      = zeros(1,length(obj.configurationSetup.activeLSNumbers));
            stimulus    = zeros(length(signal_to_play.time),length(obj.configurationSetup.activeLSNumbers));
            signal_run  = zeros(length(signal_to_play.time),length(obj.configurationSetup.activeLSNumbers));
            levels      = zeros(1,max(obj.configurationSetup.activeLSNumbers));
            
            %% VBAP balance receive the correct SPL from the two speakers
            
            angle_space = 360/length(obj.configurationSetup.ls_dir); %
            %Clockwise adjust to match odeon
            
            if obj.configurationSetup.angle> 90 && obj.configurationSetup.angle<= 180 %Desired presentation to sequential LS
                if  obj.configurationSetup.angle<= 135
                    %     s1 = 1; s2 = 19;
                    s1 = obj.configurationSetup.activeLSNumbers(obj.configurationSetup.ls_dir(:,1)==180);
                    s2 = obj.configurationSetup.activeLSNumbers(obj.configurationSetup.ls_dir(:,1)==90);
                else
                    s1 = obj.configurationSetup.activeLSNumbers(obj.configurationSetup.ls_dir(:,1)==90);
                    s2 = obj.configurationSetup.activeLSNumbers(obj.configurationSetup.ls_dir(:,1)==180);
                end
            elseif obj.configurationSetup.angle > 180 && obj.configurationSetup.angle<= 270
                if  obj.configurationSetup.angle<= 225
                    s1 = obj.configurationSetup.activeLSNumbers(obj.configurationSetup.ls_dir(:,1)==180);
                    s2 = obj.configurationSetup.activeLSNumbers(obj.configurationSetup.ls_dir(:,1)==270);
                else
                    s1 = obj.configurationSetup.activeLSNumbers(obj.configurationSetup.ls_dir(:,1)==180);
                    s2 = obj.configurationSetup.activeLSNumbers(obj.configurationSetup.ls_dir(:,1)==270);
                end
            elseif obj.configurationSetup.angle> 270 && obj.configurationSetup.angle<= 360
                if  obj.configurationSetup.angle<= 315
                    s1 = obj.configurationSetup.activeLSNumbers(obj.configurationSetup.ls_dir(:,1)==270);
                    s2 = obj.configurationSetup.activeLSNumbers(obj.configurationSetup.ls_dir(:,1)==0);
                else
                    s1 = obj.configurationSetup.activeLSNumbers(obj.configurationSetup.ls_dir(:,1)==0);
                    s2 = obj.configurationSetup.activeLSNumbers(obj.configurationSetup.ls_dir(:,1)==270);
                end
            elseif (obj.configurationSetup.angle>= 0 && obj.configurationSetup.angle<= 90)
                if  obj.configurationSetup.angle<= 45
                    s1 = obj.configurationSetup.activeLSNumbers(obj.configurationSetup.ls_dir(:,1)==0);
                    s2 = obj.configurationSetup.activeLSNumbers(obj.configurationSetup.ls_dir(:,1)==90);
                else
                    s1 = obj.configurationSetup.activeLSNumbers(obj.configurationSetup.ls_dir(:,1)==90);
                    s2 = obj.configurationSetup.activeLSNumbers(obj.configurationSetup.ls_dir(:,1)==0);
                end
            end
            
            %determine panning levels per channel
            %Coherent sum (6 dB) , so ^2 to calculate the correct energy level factor %Otherwise just sin/cos to pan law
            s1_level = (sin(abs(rem(obj.configurationSetup.angle,angle_space)/angle_space*(pi/2))))^2;
            s2_level = (cos(abs(rem(obj.configurationSetup.angle,angle_space)/angle_space*(pi/2))))^2;
            
            [s1_level, s2_level] = deal(max([s1_level s2_level]),min([s1_level  s2_level]));
            if obj.configurationSetup.level ~= 'n'
                for idx = 1:length(obj.configurationSetup.level)
                    %         [s1(idx), s2(idx) s1_level(idx) s2_level(idx)] = equal_power_pan(angle(idx));
                    if isinf(20*log10(10^((obj.configurationSetup.level(idx)/20))*s1_level(idx)))
                        levels(s1(idx)) = 0;
                    else
                        levels(s1(idx)) = 20*log10(10^((obj.configurationSetup.level(idx)/20))*s1_level(idx));
                    end
                    if isinf(20*log10(10^((obj.configurationSetup.level(idx)/20))*s2_level(idx)))
                        levels(s2(idx)) = 0;
                    else
                        levels(s2(idx)) = 20*log10(10^((obj.configurationSetup.level(idx)/20))*s2_level(idx));
                    end
                end
            end
            signal_run_ita_FILTER       = itaAudio(zeros(length(signal_to_play.time),length(obj.configurationSetup.activeLSNumbers)),signal.samplingRate,'time');
            signal_run_ita_FILTER_LEVEL = itaAudio(zeros(length(signal_to_play.time),length(obj.configurationSetup.activeLSNumbers)),signal.samplingRate,'time');
            
            for idx = 1:length(obj.configurationSetup.activeLSNumbers)
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

                signal_run(:,idx) = stimulus(:,idx).*repmat(10.^((levels(obj.configurationSetup.activeLSNumbers(idx)) - obj.configurationSetup.lsdBperVolt )./20),length(new_page(:,idx)),1);
                signal_run_ita_dB = itaAudio(signal_run,signal.samplingRate,'time');
               
                if obj.configurationSetup.default 
                    signal_run_ita_FILTER.time(:,idx) = signal_run_ita_dB.time(:,idx);
                else
                    filtered = ita_multiply_spk(signal_run_ita_dB.ch(idx),frequencyFilter.ch(obj.configurationSetup.activeLSNumbers(idx)));
                    signal_run_ita_FILTER.time(:,idx) = filtered.time(:,idx);
                end
                signal_run_ita_FILTER_LEVEL.time(:,idx) = signal_run_ita_FILTER.ch(idx).time.*(obj.configurationSetup.Level_Factor(obj.configurationSetup.activeLSNumbers(idx)));
            end

            VBAP_DS =   signal_run_ita_FILTER_LEVEL*max(DSER.time);

            for i = 1:length(obj.configurationSetup.activeLSNumbers)
                obj.VBAP_DSER.time(:,obj.configurationSetup.activeLSNumbers(i)) = VBAP_DS.time(:,i);
            end
            
            %% Combine DS ER and LR
            final_audio = ita_add(obj.VBAP_DSER,obj.Ambisonics_LR); 
            
            if length(obj.configurationSetup.LSArray) >final_audio.dimensions
                final_audio.time(:,final_audio.dimensions+1:length(obj.configurationSetup.LSArray)) = zeros(final_audio.nSamples,length(obj.configurationSetup.LSArray)-(final_audio.dimensions));
            end

            processed_audio = final_audio;

            setup_audio.iFS = obj.configurationSetup.iFS;
            setup_audio.level = obj.configurationSetup.level;
            setup_audio.angle = obj.configurationSetup.angle;
            setup_audio.signal = obj.configurationSetup.signal;
            setup_audio.VBAP_DSER = obj.VBAP_DSER;
            setup_audio.Ambisonics_LR = obj.Ambisonics_LR;

        end

        function [pansig,padsig,sig] = ring_VBAP(obj,iAngle,configurationSetup)
            %% VBAP Archontis
            % initial parameters (Here we can even move the sound source... I will
            % update this to have as an option, now is static) .Sergio.
            blocksize = obj.configurationSetup.iFS/18.3750; % (~50msec)
            if fs == 48000
                blocksize = fs/20; % (~50msec)
            elseif fs == 96000
                blocksize = fs/40; % (~50msec)
            end
            hopsize = blocksize/2; % panning hopsize (update the panning values twice per bocksize)
            ls_num = length(configurationSetup.ls_dir);
            sig = obj.configurationSetup.signal.time;
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
            
            configurationSetup.default = 1;
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
            
            configurationSetup.iFS = 44100;
            % To be used in 4 LS Hybrid mode we need to specify the LS index numbers
            configurationSetup.ls_dir = [configurationSetup.Channel_1 configurationSetup.Channel_7 configurationSetup.Channel_13 configurationSetup.Channel_19; 0 0 0 0]';
            %Right LS is 90 Degrees / Left is 270
            configurationSetup.LSArray = [mod(configurationSetup.Channel_1:configurationSetup.spaceBetweenLS:360,360) configurationSetup.spaceBetweenLS:configurationSetup.spaceBetweenLS:(configurationSetup.Channel_1-configurationSetup.spaceBetweenLS)]'; %Clockwise configuration
            %Right LS is 270 Degrees / Left is 90 counter = clockwise configuration
            % configurationSetup.LSArray = [mod(Channel_1:-spaceBetweenLS:0,360) (Channel_1-spaceBetweenLS):-spaceBetweenLS:spaceBetweenLS]';
        end
    end
end

        

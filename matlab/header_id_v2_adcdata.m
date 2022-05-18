clc 
clear
% fileID_channel0=fopen('Adc1-P-new1.bin'); %Open saved data from ADC0 quadrature 
fileID_channel0=fopen('Adc1-11-5copy.bin'); %Open saved data from ADC0 quadrature 

Data_channel0=fread(fileID_channel0,'bit8');% data in 16-bit signed format
L=length(Data_channel0); %length of recorded data
fclose(fileID_channel0);

% fileID_channel1=fopen('Ad c0-X-new1.bin'); %Open saved data from ADC0 quadrature 
fileID_channel1=fopen('Adc0-11-5copy.bin'); %Open saved data from ADC0 quadrature 
Data_channel1=fread(fileID_channel1,'bit8');% data in 16-bit signed format
L1=length(Data_channel1); %length of recorded data
fclose(fileID_channel1);
amplitude_x=Data_channel1(1:3385190);
amplitude_p=Data_channel0(1:3385190);

[phase_ori,amp_ori]=cart2pol(Data_channel1(1:3385190),Data_channel0(1:3385190));
phase_ori=wrapTo2Pi(phase_ori);

%phase and amplitude files at alice
fileID = fopen('phase_gaussian_header_10000.txt','r'); %original phase file sent from Alice
formatSpec = '%f';
alice_phase= fscanf(fileID,formatSpec);
alice_phase=(alice_phase+1)*pi;
fileID = fopen('pulse_gaussian_header_10000.txt','r'); %original phase file sent from Alice
formatSpec = '%f';
alice_amplitude= fscanf(fileID,formatSpec);
alice_amplitude=alice_amplitude+1;
first_alice=100012;
spacing_alice=50;

%extracting phase and amplitude values at alice 
for j=1:10000
    phase_alice(j)=alice_phase(first_alice+(spacing_alice*(j-1)));
    amplitude_alice(j)=alice_amplitude(first_alice+(spacing_alice*(j-1)));
    
end 

%conversion of amplitude and phase into x and p quadratures 
[alice_x, alice_p]=pol2cart(phase_alice,amplitude_alice); %original x and p sent from alice
alice_x_copy=alice_x;
alice_p_copy=alice_p;

%% header detection automation 

%parameters for header detection
k=find(amp_ori(:,1)>119);
first_ref=k(1);%find the first reference pulse in the collected data, so everything can be skipped until the following header and data block 
ref_index=first_ref;
ref_spacing=160;
sig_spacing=160;
header_yes=0;

first_ref_yes=0;
counter_ref=1;
while first_ref_yes ==0
  if(amp_ori(counter_ref)>119 && counter_ref>100)
      first_ref_yes=1;
  end
  counter_ref=counter_ref+1;
end

ref_index=counter_ref;
%locating the start of the header 
while header_yes==0
   ref_value=amp_ori(ref_index);
   if ref_value>118
       header_yes=0;
       ref_index=ref_index+ref_spacing;
   elseif amp_ori(ref_index)>98 && amp_ori(ref_index)<116 && amp_ori(ref_index+ref_spacing)>98 && amp_ori(ref_index+ref_spacing)<116
       header_yes=1;
   end
   
end 

header_index=ref_index;
ref_yes=0;

%locating the enf of the header and start of the data block
while ref_yes==0
    header_value=amp_ori(header_index);
    if header_value>98 && header_value<116
        ref_yes=0;
        header_index=header_index+ref_spacing;
    elseif amp_ori(header_index)>119 && amp_ori(header_index+ref_spacing)>119
        ref_yes=1;
    end
    
end 

counter_first_ref=1;

for counter_first_ref=1:10
    if(amp_ori(header_index-counter_first_ref)>119)
        new_header_index=header_index-counter_first_ref;
    end
end 

start_sig=new_header_index-ref_spacing/2;
start_ref=new_header_index;

%separation of x and p values at bob and shot noise values 
for j=1:10000
    bob_x_sig(j)=amplitude_x(start_sig+(j-1)*sig_spacing);
    bob_p_sig(j)=amplitude_p(start_sig+(j-1)*sig_spacing);
    bob_amp_sig(j)=amp_ori(start_sig+(j-1)*sig_spacing);
    bob_phase_sig(j)=phase_ori(start_sig+(j-1)*sig_spacing);
    
    bob_x_ref(j)=amplitude_x(start_ref+(j-1)*ref_spacing);
    bob_p_ref(j)=amplitude_p(start_ref+(j-1)*ref_spacing);
    bob_amp_ref(j)=amp_ori(start_ref+(j-1)*ref_spacing);
    bob_phase_ref(j)=phase_ori(start_ref+(j-1)*ref_spacing);
    
    shot_x(j)=amplitude_x(start_sig+40+(j-1)*sig_spacing);
    shot_p(j)=amplitude_p(start_sig+40+(j-1)*sig_spacing);
end

%% phase correction at bob 
%phase correction from reference phases 
corrected_phase=wrapTo2Pi(bob_phase_sig-bob_phase_ref);
[corrx,corrp]=pol2cart(corrected_phase, bob_amp_sig);

corrx_copy=corrx;
corrp_copy=corrp;

%removal of data values with large phase errors 
for k=1:10000

    if ((corrected_phase(k)-alice_phase(k))>=1 ||(corrected_phase(k)-alice_phase(k))<=-1 )
        corrx_copy(k)=-1000;
        corrp_copy(k)=-1000;
        alice_x_copy(k)=-1000;
        alice_p_copy(k)=-1000;   
    end
end 

corrx_final=corrx_copy(corrx_copy~=-10000);
corrp_final=corrp_copy(corrp_copy~=-10000);
alice_x_final=alice_x_copy(alice_x_copy~=-10000);
alice_p_final=alice_p_copy(alice_p_copy~=-10000);

%% excess noise and tranmittance monitoring at bob 

product_sum=sum(alice_x.*corrx)/10000;
sum_alice=sum(corrx)/10000;
sum_bob=sum(alice_x)/10000;
cov_alice_bob=product_sum-sum_alice*sum_bob; %covariance between alice and bob x quadrature values 
Va=2.12;
V_alice=var(alice_x); %alice variance calculation 
det_eff=0.5;
% shot_noise=var([shot_p shot_x]);
shot_noise=var(shot_x);

T=cov_alice_bob^2/(Va^2*det_eff); %transmittance
T_snu=T/(shot_noise);
T_1=cov_alice_bob^2/(V_alice^2*det_eff); %transmittance using variance calculated with file values 
excess_noise=(var(corrx)-det_eff*T_1*V_alice-var(shot_x)-0.1*var(shot_x))/(det_eff*T_1);
excess_noise_snu=abs(excess_noise)/sqrt(shot_noise); %excess noise in shot units 

sum_product_x=sum(alice_x.*corrx)/10000;
product_sum_x=sum(corrx.*corrx)/10000;
covab_x=sum_product_x/10000;

t_x=sqrt(2)*sum_product_x/product_sum_x;

excess_noise_x=(var(corrx)-det_eff*t_x*V_alice-sqrt(var(shot_x)));

for z=1:10000
    phase_error(z)=corrected_phase(z)-phase_alice(z);
    if(phase_error(z)>pi)
        phase_error(z)=phase_error(z)-2*pi;
    end
end 


clear
timestep= 0.01;
num_steps= 30/timestep;

   f_r_shift=3;
   f_slope= 5;
   g_r_shift=3;
   g_slope=10;
   CI_speedup=2;
   
   weightstrength_aff=8.0;
   backaff_factor = 0.5;
   
   synapse_prob=1;

   weight_strength_rec=11.0;
   
number_patterns=2; %for figure 2.B2 in paper


% CI = copetitive inhibition
CI_size=1;

% PI = Pausing the sequence inhibition
PI_size=1;


pattern_size= 2;


%legends and plotting
legend_strings= cell(number_patterns +2 ,1);
color_strings=  cell(number_patterns +2 ,1);




%weights

w_pattern_CI= ones(  CI_size,pattern_size)*2.2/pattern_size;
w_CI_pattern= ones(  pattern_size,CI_size)*(-11.0);

w_PI_pattern= ones(  pattern_size,PI_size)*(-6.0);



w_rec=binornd(1, synapse_prob, [pattern_size, pattern_size, number_patterns])*weight_strength_rec/(pattern_size*synapse_prob);

w_aff=  binornd(1, synapse_prob, [pattern_size, pattern_size, number_patterns-1])*weightstrength_aff/(pattern_size*synapse_prob);

w_aff_back=  binornd(1, synapse_prob, [pattern_size, pattern_size, number_patterns])*weightstrength_aff*backaff_factor/(pattern_size*synapse_prob);

%initial States
state0_patterns= zeros(pattern_size, number_patterns);
state0_patterns(:,1)= ones(pattern_size,1);

CI0=zeros(CI_size,1);
PI0=zeros(PI_size,1);


%timestep and storage
pattern_states= zeros(pattern_size, num_steps, number_patterns);
pattern_states(:, 1, :)= state0_patterns;


CI= zeros(CI_size, num_steps); CI(:,  1)=CI0;
PI= zeros(PI_size, num_steps); PI(:,  1)=PI0;

a=[1:1000, 1800:3000];%, 2125:3000]; %fig 2.B2

PI(:,  a)= ones(PI_size, size(a,2));


%integraion
for t=2:num_steps
   % display(t);
   %compute new state according to activation function
   %sigmoid function f paramaters
   i=1;
          input= w_rec(:,:,i)* pattern_states(:, t-1, i)+ w_aff_back(:,:,i+1)* pattern_states(:, t-1, i+1)  + w_CI_pattern * CI(:,t-1)+ w_PI_pattern * PI(:,t-1);
          new_state= 1./(1.+ exp(-(input-f_r_shift)*f_slope));
          
          pattern_states(:,t,i)= pattern_states(:,t-1,i) + timestep * (new_state - pattern_states(:,t-1,i));
   
   for i= 2: (number_patterns-1)
          input= w_rec(:,:,i)* pattern_states(:, t-1, i) + w_aff(:,:,i-1)* pattern_states(:, t-1, i-1) + w_aff_back(:,:,i+1)* pattern_states(:, t-1, i+1) + w_CI_pattern * CI(:,t-1)+ w_PI_pattern * PI(:,t-1);
          new_state= 1./(1.+ exp(-(input-f_r_shift)*f_slope));
          
          pattern_states(:,t,i)= pattern_states(:,t-1,i) + timestep * (new_state - pattern_states(:,t-1,i));
   end
   i=number_patterns;
          input= w_rec(:,:,i)* pattern_states(:, t-1, i)+ w_aff(:,:,i-1)* pattern_states(:, t-1, i-1)  + w_CI_pattern * CI(:,t-1)+ w_PI_pattern * PI(:,t-1);
          new_state= 1./(1.+ exp(-(input-f_r_shift)*f_slope));
          
          pattern_states(:,t,i)= pattern_states(:,t-1,i) + timestep * (new_state - pattern_states(:,t-1,i));
   
      
   
   sum_pattern_states= sum(pattern_states,3);
   input_CI= w_pattern_CI* (sum_pattern_states(:,t-1));
   new_state_CI= 1./(1.+ exp(-(input_CI-g_r_shift)*g_slope));
   
   CI(:,t)= CI(:,t-1) + CI_speedup*timestep*(new_state_CI - CI(:,t-1));
  
   
    
    
    
end
for i=1:number_patterns
%plotting
    %new_legend_string=['B_', num2str(i)];
    %legend_strings{i}=new_legend_string;
    if mod(i,6)== 1
        new_legend_string=['B_i'];
        legend_strings{i}=new_legend_string;
        new_color_string=['r'];
        color_strings{i}=new_color_string;
    elseif mod(i,6)== 2
        new_legend_string=['B_{i+1}'];
        legend_strings{i}=new_legend_string;
        new_color_string=['m'];
        color_strings{i}=new_color_string;
    elseif mod(i,6)== 3
        new_color_string=['y'];
        color_strings{i}=new_color_string;
    elseif mod(i,6)== 4
        new_color_string=['y'];
        color_strings{i}=new_color_string;
    elseif mod(i,6)== 5
        new_color_string=['m'];
        color_strings{i}=new_color_string;
    else
        new_color_string=['k'];
        color_strings{i}=new_color_string;
    end
end
color_strings{number_patterns+1}= ['c'];
color_strings{number_patterns+2}= ['b'];
legend_strings{number_patterns+1}= ['I_C'];
legend_strings{number_patterns+2}= ['I_G'];

X= timestep:timestep:timestep*(num_steps);
Y=zeros(num_steps, number_patterns + 2);
Y(:,1:number_patterns)= [mean(pattern_states)];
Y(:,number_patterns+1)= CI;
Y(:,number_patterns+2)= PI;
%Y= [A;B;C;D;E;CI;PI];


    
    

%%%% Figure for 3 patterns -> 2.B2
hfig = figure(1)
set(hfig, 'Position', [0 0 100 200])
h(1) = subplot(2,1,1)
H1 = plot(X, Y(:,1:number_patterns+1));
AX=gca;
    %h_legend=legend(legend_strings);
    h_legend=legend(legend_strings(1:number_patterns+1));
    set(H1, {'color'},color_strings(1:number_patterns+1));
    %set(H1, {'marker'},{'*'; '*'; '*'; '*'; '*'; '*'; '+'; 'o'})
    %set(H1, {'markers'}, {10; 10; 10; 10; 10; 10 ; 10; 10})
    set(H1, {'LineStyle'},{'-'; '--'; '-.'})
    set(H1, {'LineWidth'},{6; 6; 6})
    %title('Activity evolution', 'Interpreter', 'latex','FontSize', 45)
    set(h_legend,'Interpreter', 'latex','FontSize',50);
    legend('Location', 'NorthEast');
    %, 'Mmax / (p_aff+p_rec)'
    %xlabel('time','Interpreter', 'latex', 'FontSize', 40);
    %xlim([10^3, 10^17]);
    ylabel( '\textbf{Rate}','Interpreter', 'latex', 'FontSize', 70);
    xlim([0 30])
    set(AX, 'FontSize', 50)   
    set(AX,'XTickLabel',[]);
    set(AX,'YTickLabel',[]);

    set(h(1),'xticklabel',[]);
    SP=10; %your point goes here
    line([SP SP],get(AX,'YLim'),'Color',[0 0 1],'LineWidth', 6)
    SP=18; %your point goes here
    line([SP SP],get(AX,'YLim'),'Color',[0 0 1],'LineWidth', 6)
        
    x = [0.34 0.34];
    y = [0.91 0.955];
    annotation('textarrow',x,y,'LineWidth', 5,'String','(i)','FontSize', 50);
    
    x = [0.368 0.388];
    y = [0.43 0.46];
    annotation('textarrow',x,y,'LineWidth', 5, 'String','(ii)','FontSize', 50);
    
    x = [0.43 0.4];
    y = [0.64 0.67];
    annotation('textarrow',x,y,'LineWidth', 5,'String','(iii)','FontSize', 50);
    
    x = [0.5 0.5];
    y = [0.89 0.94];
    annotation('textarrow',x,y,'LineWidth', 5,'String','(iv)','FontSize', 50);
    
    x = [0.615 0.595];
    y = [0.43 0.46];
    annotation('textarrow',x,y,'LineWidth', 5,'String','(v)','FontSize', 50);
    
    x = [0.62 0.62];
    y = [0.91 0.955];
    annotation('textarrow',x,y,'LineWidth', 5,'String','(vi)','FontSize', 50);

h(2) = subplot(2,1,2)
H1 = plot(X, Y(:,number_patterns+2:number_patterns+2));
 AX=gca;
    %h_legend=legend(legend_strings);
    h_legend=legend(legend_strings(number_patterns+2:number_patterns+2));
    set(H1, {'color'},color_strings(number_patterns+2:number_patterns+2));
    %set(H1, {'marker'},{'*'; '*'; '*'; '*'; '*'; '*'; '+'; 'o'})
    %set(H1, {'markers'}, {10; 10; 10; 10; 10; 10 ; 10; 10})
    set(H1, {'LineStyle'},{'--'})
    set(H1, {'LineWidth'},{6})
    %title('Ensemble activity evolution', 'Interpreter', 'latex','FontSize', 50)
    set(h_legend,'Interpreter', 'latex','FontSize',50);
    legend('Location', 'NorthEast');
    %, 'Mmax / (p_aff+p_rec)'
    xlim([0 60])
    xlabel('\textbf{Time}','Interpreter', 'latex', 'FontSize', 70);
    %xlim([10^3, 10^17]);
    ylabel( '\textbf{Rate}','Interpreter', 'latex', 'FontSize', 70);
    set(AX, 'FontSize', 50)
    set(AX,'XTickLabel',[]);
    set(AX,'YTickLabel',[]);
    linkaxes(h,'x');
    
    
    
pos=get(h,'position');
bottom=pos{2}(2);
top=pos{1}(2)+pos{1}(4);
plotspace=top-bottom;
pos{2}(4)=plotspace/2-0.04;
pos{1}(4)=plotspace/2+0.04;
pos{1}(2)=bottom+plotspace/2;

set(h(1),'position',pos{1});
set(h(2),'position',pos{2});
set(h(1),'ytick',[-0.5 0 0.5 1]);
%set(h(2),'YAxisLocation','right')



    

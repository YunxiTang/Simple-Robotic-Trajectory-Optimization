classdef sol_mdl < handle
    %SOL_MDL used to store all the iteration results in cell form
    
    properties
        num_trials = 0,
        XTrials,
        UTrials,
        time_bar,
        Defects
    end
    
    methods
        function obj = sol_mdl(iter_max)
            %SOL_MDL Construct an instance of this class
            %   Detailed explanation goes here
            disp('[INFO]: Creating a sol_mdl.');
            obj.XTrials = cell(iter_max,1);
            obj.UTrials = cell(iter_max,1);
            obj.Defects = cell(iter_max,1);
        end
        
        function [] = store_trial(obj,Xtrial,Utrial)
            % store trails
            obj.num_trials = obj.num_trials + 1;
            i = obj.num_trials;
            obj.XTrials{i} = Xtrial; 
            obj.UTrials{i} = Utrial;
        end
        
        function [] = store_time(obj, time_seq)
            % store time sequence
            obj.time_bar = time_seq;
        end
        
        function [x_dft] = store_defect(obj,trail)
            % store state defect at segment nodes
            % extract from trails
            [S, ~ , Nx, ~] = size(trail);
            x_dft = zeros(S-1, Nx, 1);
            for i=1:S-1
                x_dft(i,:,:) = trail(i,end,:) - trail(i+1,1,:);
            end
            obj.Defects{obj.num_trials} = x_dft;
        end
        
        function [flatten_traj] = flatted_traj(obj, trail)
            % flatted state trajectory
            [S,LL,Nx] = size(trail);
            L = LL - 1;
            flatten_traj = zeros(S*L+1,Nx,1);
            k = 1;
            for i=1:S
                if i==S
                    EN_L = L + 1;
                else
                    EN_L = L;
                end
                for j=1:EN_L
                    flatten_traj(k,:,:) = trail(i,j,:);
                    k = k + 1;
                end
            end
        end
        
        function [flatten_u] = flatted_u(obj, u_trail)
            % flatted control trajectory
            [S,L,Nu] = size(u_trail);
            flatten_u = zeros(S*L,Nu,1);
            k = 1;
            for i=1:S 
                for j=1:L
                    flatten_u(k,:,:) = u_trail(i,j,:);
                    k = k + 1;
                end
            end
        end
    end
end


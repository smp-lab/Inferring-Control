classdef VestibularTFs_DS
    %VestibularTFs: Vestibular transfer functions object 
        % Includes the transfer function
        % parameters and plotting for the Canal and Otolith afferents
        % (regular and irregular) responses to motion and EVS inputs.
        %
        % Note: Motion inputs should be transformed into the specific
        % canal/otolith axes prior to this. To do so, use:
        %     canalObj, and
        %     otolithObj


    properties   
        % Vestibular Motion-to-Afferent transfer function parameter values
        % from Schneider et al. (2015)
        % doi: 10.1523/JNEUROSCI.3841-14.2015
        %{
          Motion to Canal Afferent transfer function take the form: 
            
                Hcanal(s) = K *   s*(s+1/T1)
                               ----------------
                              (s+1/Tc) * (s+1/T2)

          Motion to Otolith Afferent Transfer Function Take the Form: 
            
                Hoto(s) = K *   s^k1 * (1 + a*s)^k2
                                 ----------------------
                                       (1 + b*s)

        %}

        Canal_Motion_Regular = [0.0175 0.0027 5.7 2.83]; % TF parameter values [T1 T2 Tc K]
        Canal_Motion_Irregular = [0.03 0.0006 5.7 27.09]; % TF parameter values [T1 T2 Tc K]
        Otolith_Motion_Regular = [0.0643 2.208 0.0138 0.0255 59.0106] %  TF parameter values [k1 k2 a b K]
        Otolith_Motion_Irregular = [0.3084 2.6834 0.0136 0.0318 112.7417] % TF parameter values [k1 k2 a b K]
        
        Otolith_BP = [0.62 0.016 22.05 1.53] %[5.33 0.66 13.2 0.4] % [TL Ts Ta K] REMOVE???

        % Vestibular EVS-to-Afferent transfer function parameters
        % from Kwan et al. (2019)
        % doi:10.1038/s41467-019-09738-1
        %{
         
            EVS to Canal Afferent transfer function take the form:

                                       K(s+b1)(s+b2)
                   H_Canal_EVS  =   -------------------
                                       (s+a1)(s+a2)

            EVS to Otolith Afferent transfer function take the form:

                                       s^k1(1 + b*s)^k2
               H_EVS_regular(s) =  K ---------------------
                                           (1 + a*s)
        %}

        Canal_EVS_Regular = [26 188 47 2578 98]; % TF parameter values [b1 b2 a1 a2 K]
        Canal_EVS_Irregular = [12 136 18 2739 418]; % TF parameter values [b1 b2 a1 a2 K]
        Otolith_EVS_Regular = [0.09 1.30 0.019 0.014 3.15]; % TF parameter values [k1 k2 b a K]
        Otolith_EVS_Irregular = [0.12 1.46 0.009 0.009 13.8]; % TF parameter values [k1 k2 b a K]

    end 

    methods
        function [H_Canal_Motion_Regular, H_Canal_Motion_Irregular, H_Canal_Motion_Mixed] = canalMotionTFs(obj)
            % Generate transfer functions for canal motion-to-afferent responses
            %{
                Semicircular Canal: Angular Velocity to Primary Canal
                Afferent Transfer Function. This assumes a 3:1 ratio of
                regular and irregular afferent population. 
                
                Canal to Afferent Transfer Function Take the Form: 
            
                Hcanal(s) = K *   s*(s+1/T1)
                               ----------------
                              (s+1/Tc) * (s+1/T2)
            
            
                Returns: 
                1) Transfer Function for Regular (scaled) 
                2) Transfer Function for Irregular (scaled)
                3) Transfer Function for Mixed (regular + irregular)
            %}
            
            Regular_T1 = obj.Canal_Motion_Regular(1); 
            Regular_T2 = obj.Canal_Motion_Regular(2); 
            Regular_Tc = obj.Canal_Motion_Regular(3); 
            Regular_K = obj.Canal_Motion_Regular(4); 
            
            Irregular_T1 = obj.Canal_Motion_Irregular(1);
            Irregular_T2 = obj.Canal_Motion_Irregular(2); 
            Irregular_Tc = obj.Canal_Motion_Irregular(3); 
            Irregular_K = obj.Canal_Motion_Irregular(4); 
            
            H_Canal_Motion_Regular = Regular_K * tf( [1, 1./Regular_T1, 0], [1, (1./Regular_T2)+(1./Regular_Tc), (1./Regular_Tc)*(1./Regular_T2)] );
            H_Canal_Motion_Irregular = Irregular_K * tf( [1, 1./Irregular_T1, 0], [1, (1./Irregular_T2)+(1./Irregular_Tc), (1./Irregular_Tc)*(1./Irregular_T2)] );           
            H_Canal_Motion_Mixed = (H_Canal_Motion_Regular * 0.25 * 3) + (H_Canal_Motion_Irregular* 0.25);
        end
        
       
        
        function [Hoto_Regular, Hoto_Irregular, Hoto_Mixed] = otolithMotionTFs(obj) 
            % Generate transfer functions for otolith motion-to-afferent responses
             %{
                Otolith: Acceleration (in G) to Primary Otolith
                Afferent Transfer Functions.  
                
                * NOTE Fractional Order Transfer Function Requires FOMCON
                toolbox and fotf library  
            
                Otolith to Afferent Transfer Function Take the Form: 
            
                Hoto(s) = K *   s^k1 * (1 + a*s)^k2
                                 ----------------------
                                       (1 + b*s)
            
                Returns: 
                1) Transfer Function for Regular (scaled) 
                2) Transfer Function for Irregular (scaled)
                3) Transfer Function for Mixed (regular + irregular)
            %}
           
            
            Hoto_Regular = obj.otoTF(obj.Otolith_Motion_Regular);
            Hoto_Irregular = obj.otoTF(obj.Otolith_Motion_Irregular);

            Hoto_Mixed = (0.25 * 3 * Hoto_Regular) + (0.25 * Hoto_Irregular);
        end 
        
    end

  methods (Static) 
        function G = otoTF(parameters)
            % internal method to generate fractional order transfer function
            k1 = parameters(1);
            k2 = parameters(2);
            a = parameters(3);
            b = parameters(4);
            K = parameters(5); 
            
            syms s
            f = ((1 + a*s)^(k2));
            T = vpa(taylor(f, s, 'order', 8, 'OrderMode','absolute', 'ExpansionPoint',0), 10);
            f2 = s^k1;
            
            str2 = char(vpa(1 + b*s,10));
            
            T_str = char(T);
            n = strfind(T_str,'^');
            for i=1:length(n)
                T_str = [T_str(1:n(i)) '{' T_str(n(i)+1) '}' T_str(n(i)+2:end)];
                n = n + 2;
            end
            
            str1 = regexprep( char(f2) ,'[(]','{' );
            str1 = regexprep( str1 ,'[)]','}' );
            
            G = K * (fotf(char(T), '1') * fotf(str1, '1') * fotf('1', str2));
        end
  end
end
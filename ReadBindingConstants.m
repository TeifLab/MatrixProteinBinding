function [TF_profile] = ReadBindingConstants(varargin)

% function uses read in DNA_seq to calculate binding constants for a transcription factor
% using PWM method:
%
% Default is one argument in, and PWM matrix for CTCF
%
% call PWM method for base-pair level binding constants

global K0


% set some defaults

g = 1;
matrix_file = 'CTCF_matrix_Orlov_txt';

if (nargin>1)
   g = varargin{3};
   matrix_file = varargin{2};
end

DNA_seq = varargin{1};

DNA_seq = abs(DNA_seq);

[K_pwm,K_rev,~,~] = PWM_trap(DNA_seq,K0(g),1.5,matrix_file);

TF_profile = K_pwm + K_rev;

end

function doPitchShift(shift_st)
% shift_st: the amount  of shift in semitones (0: noshift)
% Original path: d:\speechres\kape\ADAPT_VC_RATIO\BIN\release\TransShiftMex.mexw32
if ~exist('shift_st')
    shift_st = 3;
end

TSM_DAF_PATH = '../Audapter-DAF/BIN/release/';

% origPath = which('TransShiftMex');

addpath(TSM_DAF_PATH);
TransShiftMex_PS(3, 'srate', 16000);
TransShiftMex_PS(3, 'framelen', 32);
TransShiftMex_PS(3, 'bpitchshift', 1);
TransShiftMex_PS(3, 'pitchshiftratio', 2 ^ (shift_st / 12));
TransShiftMex_PS(3, 'pvocframelen', 256);
TransShiftMex_PS(3, 'pvochop', 64);

TransShiftMex_PS(1);
% pause(10);
% TransShiftMex_PS(2);

% addpath(origPath);
return
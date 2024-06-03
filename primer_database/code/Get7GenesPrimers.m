for i = 1:length(FGFR_rev)
   FGFR_rev_OFP{i} = FindFP(FGFR_rev{i},[11.5 12.5],60); 
   FGFR_rev_OFP_Seq{i} = FGFR_rev_OFP{i}{1};
   FGFR_rev_OFP_Eny{i} = FGFR_rev_OFP{i}{3};
   FGFR_rev_OFP_Pos{i} = FGFR_rev_OFP{i}{2};
   FGFR_rev_Seq_Len{i} = length(FGFR_rev{i})-FGFR_rev_OFP_Pos{i}(2);
end

% FGFR_rev

FGFR_rev_OFP_Seq = FGFR_rev_OFP_Seq';
FGFR_rev_OFP_Eny = FGFR_rev_OFP_Eny';
FGFR_rev_OFP_Pos= FGFR_rev_OFP_Pos';
FGFR_rev_Seq_Len = FGFR_rev_Seq_Len';
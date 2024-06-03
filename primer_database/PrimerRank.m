function [RankValue] = PrimerRank(dG,GCcontent)

if (dG >= -12.5) & (dG <= -10.5)
    RankValue = 1000;
else
    RankValue = -abs(dG+11.5);
end

if (GCcontent >= 0.4) & (GCcontent <= 0.6)
    RankValue = RankValue + 100 - abs(GCcontent-0.5)*100;
else
    RankValue = RankValue - abs(GCcontent-0.5)*100;
end

% 
fn = 'dataProcessed/dataset/dataNeurophys_SAT.mat';
dNphys = load(fn);

unitStats = dNphys.unitStats;

% Col names, array vars
colNamesArrayVars =  {
    {'Baseline_SAT_Effect'           }    {'Signed logical'                                    }
    {'Baseline_Stats'                }    {'Acc_Mean  Fast_Mean  Acc_SD  Fast_SD'              }
    {'BuildupActivity'               }    {'Acc_Correct  Fast_Correct  Acc_Error  Fast_Error'  }
    {'ChoiceErrorSignal_CheckRatioDT'}    {'Ratio  tSignal_Fast_ShortISI  tSignal_Fast_LongISI'}
    {'ChoiceErrorSignal_Magnitude'   }    {'Accurate  Fast'                                    }
    {'ChoiceErrorSignal_Time'        }    {'t_Acc_Start  t_Acc_End  t_Fast_Start  t_Fast_End'  }
    {'NormalizationFactor'           }    {'All  Vis  Move  Reward'                            }
    {'TargetSelectionTime'           }    {'Accurate  Fast'                                    }
    {'TimingErrorSignal_Magnitude'   }    {'Accurate  Fast'                                    }
    {'TimingErrorSignal_Time'        }    {'t_Acc_Start  t_Acc_End  t_Fast_Start  t_Fast_End'  }
    {'VisualResponse_Latency'        }    {'Time from array'                                   }
    {'VisualResponse_Magnitude'      }    {'Accurate  Fast'                                    }
    {'VisualResponse_SAT_Effect'     }    {'Logical'                                           }
};
% Map col name arrays to col names and descriptiomns
colNamesVarMap = containers.Map();
colNamesVarMap('Baseline_SAT_Effect') = ...
    {{'BaselineSatEffect'} 
     {'Baseline_SAT_Effect Signed logical'}};
colNamesVarMap('Baseline_Stats') = ...
   {{'BaselineAccMean','BaselineFastMean','BaselineAccStd','BaselineFastStd'}
    {'Baseline_Stats Acc_Mean','Baseline_Stats Fast_Mean','Baseline_Stats Acc_SD','Baseline_Stats Fast_SD'}};
colNamesVarMap('BuildupActivity') = ...
   {{'BuildupActivityAccCorrect','BuildupActivityFastCorrect','BuildupActivityAccError','BuildupActivityFastError'} 
    {'BuildupActivity Acc_Correct','BuildupActivity Fast_Correct','BuildupActivity Acc_Error','BuildupActivity Fast_Error'}};
colNamesVarMap('ChoiceErrorSignal_CheckRatioDT') = ...
   {{'ErrorChoiceSignalDtRatio','ErrorChoiceSignalFastShortIsi','ErrorChoiceSignalFastLongIsi'} 
    {'ChoiceErrorSignal_CheckRatioDT Ratio','ChoiceErrorSignal_CheckRatioDT tSignal_Fast_ShortISI','ChoiceErrorSignal_CheckRatioDT tSignal_Fast_LongISI'}};
colNamesVarMap('ChoiceErrorSignal_Magnitude') = ...
   {{'ErrorChoiceSignalMagAcc','ErrorChoiceSignalMagFast'} 
    {'ChoiceErrorSignal_Magnitude Accurate','ChoiceErrorSignal_Magnitude Fast'}};
colNamesVarMap('ChoiceErrorSignal_Time') = ...
   {{'ErrorChoiceSignalTimAccStart','ErrorChoiceSignalTimAccEnd','ErrorChoiceSignalTimFastStart','ErrorChoiceSignalTimFastEnd'} 
    {'ChoiceErrorSignal_Time t_Acc_Start','ChoiceErrorSignal_Time t_Acc_End','ChoiceErrorSignal_Time t_Fast_Start','ChoiceErrorSignal_Time t_Fast_End'}};
colNamesVarMap('NormalizationFactor') = ...
   {{'NormalizationFactorAll','NormalizationFactorVisual','NormalizationFactorMove','NormalizationFactorReward'} 
    {'NormalizationFactor All','NormalizationFactor Vis','NormalizationFactor Mov','NormalizationFactor Reward'}};
colNamesVarMap('TargetSelectionTime') = ...
   {{'TargetSelectionTimeAcc','TargetSelectionTimeFast'} 
    {'TargetSelectionTime Accurate','TargetSelectionTime Fast'}};
colNamesVarMap('TimingErrorSignal_Magnitude') = ...
   {{'ErrorTimingSignalMagAcc','ErrorTimingSignalMagFast'} 
    {'TimingErrorSignal_Magnitude Accurate','TimingErrorSignal_Magnitude Fast'}};
colNamesVarMap('TimingErrorSignal_Time') = ...
   {{'ErrorTimingSignalTimAccStart','ErrorTimingSignalTimAccEnd','ErrorTimingSignalTimFastStart','ErrorTimingSignalTimFastEnd'} 
    {'TimingErrorSignal_Time t_Acc_Start','TimingErrorSignal_Time t_Acc_End','TimingErrorSignal_Time t_Fast_Start','TimingErrorSignal_Time t_Fast_End'}};
colNamesVarMap('VisualResponse_Latency') = ...
    {{'VisualResponseLatency'} 
     {'VisualResponse_Latency Time from array'}};
colNamesVarMap('VisualResponse_Magnitude') = ...
   {{'VisualResponseMagAcc','VisualResponseMagFast'} 
    {'VisualResponse_Magnitude Accurate','VisualResponse_Magnitude Fast'}};
colNamesVarMap('VisualResponse_SAT_Effect') = ...
    {{'VisualResponseSatEffect'} 
     {'VisualResponse_SAT_Effect Logical'}};

ustats2 = table();

colNames = unitStats.Properties.VariableNames;
for ii = 1:numel(colNames)
    colName = colNames{ii};
    vals = unitStats{:,colName};
    newColNameDesc = colNamesVarMap(colName);
    temp = table();
    if size(vals,2) == 1
        temp.tempName = vals;
        temp.Properties.VariableNames = newColNameDesc{1};
        temp.Properties.VariableDescriptions = newColNameDesc{2};
    else
        temp = array2table(vals);
        temp.Properties.VariableNames = newColNameDesc{1};
        temp.Properties.VariableDescriptions = newColNameDesc{2};
    end
    ustats2 = [ustats2 temp];
    if ii == 1
        ustats2.Properties.RowNames = unitStats.Properties.RowNames;
    end
end

dNphys.unitStatsOld = dNphys.unitStats;
dNphys.unitStats = ustats2;

save(fn,'-v7.3','-struct','dNphys');






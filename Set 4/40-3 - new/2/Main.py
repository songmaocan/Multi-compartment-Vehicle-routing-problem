"""
Date: 2021-4-17
Function: ALR based method
Version: 1.1
Update: Number of vehicles and time intervals is not need, add solution state,output GAMS file
"""
from Model import LR
import time
def main():
    start_time=time.time()
    mod=LR()
    mod.g_solve_CSMC_by_LR()
    end_time=time.time()
    spend_time=end_time-start_time
    mod.output_to_file(spend_time)
    print("CPU running time {} min".format(spend_time/60))

if __name__=="__main__":
    main()
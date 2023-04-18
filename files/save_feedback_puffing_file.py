#!/usr/bin/python

# Function definition is here

import os
from routines.globals		import DEBUG

#=========================================================
# This routine to write data to the H5DF mesh file
#=========================================================

def save_feedback_puffing_file(FeedbackFilePath, FeedbackPuffing): 

	if(DEBUG > 0): print("save_feedback_puffing_file: Saving to ",FeedbackFilePath+"feedback_data.txt")

	try:
		os.mkdir(FeedbackFilePath)
	except OSError:
		pass


#	Feedback option input file
	
	fid = open(FeedbackFilePath+"feedback_data.txt",'w')

	fid.write("{:0.4e}    ! n target\n".format(FeedbackPuffing.nTarget))
	fid.write(    "{:d}    ! i target \n".format(FeedbackPuffing.iTarget+1))
	fid.write(    "{:d}    ! j target \n".format(FeedbackPuffing.jTarget+1))
	fid.write(    "{:d}    ! k target \n".format(FeedbackPuffing.kTarget+1))
	fid.write(    "{:f}    ! Gain\n".format(FeedbackPuffing.Gain))
	fid.write("{:0.4e}    ! tau_i\n".format(FeedbackPuffing.Tau_i))
	fid.write("{:0.4e}    ! tau_d\n".format(FeedbackPuffing.Tau_d))
	fid.write("{:0.4e}    ! max puff\n".format(FeedbackPuffing.MaxPuff))
	fid.write("{:0.4e}    ! min puff\n".format(FeedbackPuffing.MinPuff))

	fid.close()

	return

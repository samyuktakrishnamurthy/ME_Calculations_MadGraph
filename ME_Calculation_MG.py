#This code reads a minitree and does the following
#calculate the FV
#Puts the FV through MadGraph
#Calculates the ME for the ???
#creates a new minitree with these new variables

import subprocess
import numpy as np
from ROOT import *
import os
import shutil
from datetime import datetime
start=datetime.now()

#Minitree we are reading from

# defining th paths for different directories
home2 = "/scratch/skrishna"
scratchArea = "/scratch/"

#os.mkdir("skrishna")
os.chdir(scratchArea)

#Defining all the paths to the calculators, based on final state
ME_ggHZZ_4l = "./ggHZZ_0jet_rclsa/SubProcesses/PV0_0_1_gg_epemepem"
ME_ggHZZ_2l2l = "./ggHZZ_0jet_rclsa/SubProcesses/PV0_0_2_gg_epemmupmum"

ME_ggnoHZZ_4l = "./ggnoHZZ_0jet_rclsa/SubProcesses/PV0_0_1_gg_epemepem"
ME_ggnoHZZ_2l2l = "./ggnoHZZ_0jet_rclsa/SubProcesses/PV0_0_2_gg_epemmupmum"

ME_qqZZ_4l = "./qqZZ_0jet_rclsa/SubProcesses/P0_uux_epemepem"
ME_qqZZ_2l2l = "./qqZZ_0jet_rclsa/SubProcesses/P2_uux_epemmupmum"

#Defining the TLorentzVectors
lep1Z1 = TLorentzVector()
lep2Z1 = TLorentzVector() 
lep1Z2 = TLorentzVector()
lep2Z2 = TLorentzVector()
jet1 = TLorentzVector() 
jet1_corr = TLorentzVector() 
gluon1 = TLorentzVector()
gluon2 = TLorentzVector()
lep2Z2_corr = TLorentzVector()

#module that calculates the ME
def Calculate_ME(event_num, loop_size, pathToFile, process):

    i = event_num
    #define home
    TarArea = scratchArea+"tmp_"+process+"_"+str(i)
    myfile = TFile(pathToFile)
    #get the minitree from the TFile and loop over the entries
    mytree = myfile.Get('nominal')
    entries = mytree.GetEntriesFast()  
    #output tree with new branches
    output_tree = process+"_"+str(i)+ ".root"
    f = TFile(output_tree,'RECREATE')
    
    #creating a new directory in the scratch area and moving there
    if (os.path.isdir("tmp_"+process+"_"+str(i))==True):
        file_copy = "tmp_"+process+"_"+str(i)
        shutil.rmtree(file_copy) 
    #make a directory and work in here for each i (basically each job)
    #Now you are in that specific area
    os.mkdir("tmp_"+process+"_"+str(i))
    os.chdir("tmp_"+process+"_"+str(i))

    #defining i/o files for check_sa.f to read and write to prevent overwriting
    ps_input = TarArea + "/input"+process+str(i)+".input"
    results = TarArea + "/output"+process+str(i)+".dat"
    
    #scping the MG code over to scratch and untaring the files
    pathToTarMG = "abc-at13:/gpfs3/umass/HZZ/MG5_aMC_v2_7_2/"

    subprocess.call(["scp", pathToTarMG+"Tar_ggHZZ_0jet_rclsa.tar.gz","./"])
    subprocess.call(["scp", pathToTarMG+"Tar_ggnoHZZ_0jet_rclsa.tar.gz","./"])
    subprocess.call(["scp", pathToTarMG+"Tar_qqZZ_0jet_rclsa.tar.gz","./"])

    subprocess.call(["tar","-xzf","Tar_ggHZZ_0jet_rclsa.tar.gz"])
    subprocess.call(["tar","-xzf","Tar_ggnoHZZ_0jet_rclsa.tar.gz"])
    subprocess.call(["tar","-xzf","Tar_qqZZ_0jet_rclsa.tar.gz"])

    #defining all the required variables for new branches
    MG_ggHZZ_ME = np.array([0.])
    MG_ggnoHZZ_ME = np.array([0.])
    MG_NLO_qqZZ_ME = np.array([0.])

    #cloning the old tree to get a new tree
    newTree = mytree.CloneTree(0)
    
    #adding all the new ME branches
    newTree.Branch("MG_ggHZZ_ME", MG_ggHZZ_ME, 'MG_ggHZZ_ME/D')
    newTree.Branch("MG_ggnoHZZ_ME", MG_ggnoHZZ_ME, 'MG_ggnoHZZ_ME/D')
    newTree.Branch("MG_NLO_qqZZ_ME", MG_NLO_qqZZ_ME, 'MG_NLO_qqZZ_ME/D')
    
    total = mytree.GetEntries()
    #print total
    
    #compile and calculate the ME
    #subprocess.call(["make","clean"])
    #subprocess.call(["make","check"])

    #Loop over entries in the trees
    i_events = event_num*loop_size
    for jentry in range(i_events, i_events+loop_size):
        ientry = mytree.LoadTree(jentry)
        if ientry < 0:
            break
        nb = mytree.GetEntry(jentry)
        if nb<=0:
            continue 
        
        #Prepare 4-vectors for the two leptons of the first Z
        lep1Z1.SetPtEtaPhiM(mytree.lepton_pt[0], mytree.lepton_eta[0], mytree.lepton_phi[0], mytree.lepton_m[0])
        lep2Z1.SetPtEtaPhiM(mytree.lepton_pt[1], mytree.lepton_eta[1], mytree.lepton_phi[1], mytree.lepton_m[1])
        #Prepare 4-vectors for the two leptons of the second Z
        lep1Z2.SetPtEtaPhiM(mytree.lepton_pt[2], mytree.lepton_eta[2],mytree.lepton_phi[2], mytree.lepton_m[2])
        lep2Z2.SetPtEtaPhiM(mytree.lepton_pt[3], mytree.lepton_eta[3],mytree.lepton_phi[3], mytree.lepton_m[3])
    
        #set Default values for ME
        MG_ggHZZ_ME_value = 0.1
        MG_ggnoHZZ_ME_value = 0.1
        MG_NLO_qqZZ_ME_value = 0.1

        ########################
        #   inclusive JET ME   #
        ########################
        
        if mytree.n_jets >=0: # This allows basically all the evnents to pass 
            
            #Reconstructing 0 jet higgs
            vecZ1= (lep1Z1 + lep2Z1)
            vecZ2= (lep1Z2 + lep2Z2)
            vecHiggs = vecZ1 + vecZ2
            #Boosting the leptons to the CM frame
            lep1Z1.Boost(-vecHiggs.BoostVector())
            lep1Z2.Boost(-vecHiggs.BoostVector())
            lep2Z1.Boost(-vecHiggs.BoostVector())
            lep2Z2.Boost(-vecHiggs.BoostVector())
            #Setting the gluons with pz as half the energy of the higgs
            gluon1.SetPxPyPzE(0,0,vecHiggs.M()/2, vecHiggs.M()/2)
            gluon2.SetPxPyPzE(0,0,-vecHiggs.M()/2, vecHiggs.M()/2)

            #Correcting the fourth lepton fourvectors to impose energy-momentum conservation
            #was causing Energy-momentum conservation issues, although the difference was of the order 10^-13 or so
            lep2Z2_corr_px = - (lep1Z1.Px()+lep2Z1.Px()+lep1Z2.Px())
            lep2Z2_corr_py = - (lep1Z1.Py()+lep2Z1.Py()+lep1Z2.Py())
            lep2Z2_corr_pz = - (lep1Z1.Pz()+lep2Z1.Pz()+lep1Z2.Pz())
            lep2Z2_corr_e = np.sqrt(lep2Z2_corr_px**2+lep2Z2_corr_py**2+lep2Z2_corr_pz**2)
            lep2Z2_corr.SetPxPyPzE(lep2Z2_corr_px, lep2Z2_corr_py, lep2Z2_corr_pz, lep2Z2_corr_e);

            #File to which the FV is written that is the input for the check_sa.f
            FV_file = open(ps_input, "w")
            #Output file for the Fourvectors
            if lep1Z1.Px()==lep1Z1.Px():
                FV_file.write('%f %f %f %f \n' %(gluon1.E(),gluon1.Px(), gluon1.Py(), gluon1.Pz()))
                FV_file.write('%f %f %f %f \n' %(gluon2.E(),gluon2.Px(), gluon2.Py(), gluon2.Pz()))           
                FV_file.write('%f %f %f %f \n' %(lep1Z1.E(),lep1Z1.Px(), lep1Z1.Py(), lep1Z1.Pz()))
                FV_file.write('%f %f %f %f \n' %(lep2Z1.E(),lep2Z1.Px(), lep2Z1.Py(), lep2Z1.Pz()))
                FV_file.write('%f %f %f %f \n' %(lep1Z2.E(),lep1Z2.Px(), lep1Z2.Py(), lep1Z2.Pz())) 
                FV_file.write('%f %f %f %f \n' %(lep2Z2.E(),lep2Z2.Px(), lep2Z2.Py(), lep2Z2.Pz()))
            FV_file.close()     
                
            ########################
            #     ggHZZ ME-inc     #
            ########################     
            os.chdir(TarArea) 
            #what type of event is it -going to the right directory based on final state
            
            if(mytree.event_type<2):
                os.chdir(ME_ggHZZ_4l)
                subprocess.call(["./check",ps_input,results])
            if(mytree.event_type>=2):
                os.chdir(ME_ggHZZ_2l2l)
                subprocess.call(["./check",ps_input,results])
            #Read the ME from the results file
            file = open(results, "r")                
            outputLines = file.readlines()
            #getting the right line from the results file, and reading the ME value, 
            imp_line = outputLines[8]
            ME_list1 = imp_line.split('         ')                
            #print type(ME_list[1])
            MG_ggHZZ_ME_value = float(ME_list1[1])
            os.remove(results)
            #going back to home directory to start another calculation
            os.chdir(TarArea)
           
            ########################
            #     ggZZ ME-inc       #
            ########################
            
            #Going to the right directory
            if(mytree.event_type<2):
                os.chdir(TarArea)
                os.chdir(ME_ggnoHZZ_4l)
                subprocess.call(["./check",ps_input,results])
            if(mytree.event_type>=2):
                os.chdir(TarArea)
                os.chdir(ME_ggnoHZZ_2l2l)
                subprocess.call(["./check",ps_input,results])
            #Read the ME from the results file
            file = open(results, "r")
            outputLines = file.readlines()
            #getting the right line from the results file, and reading the ME value,
            imp_line = outputLines[8]
            ME_list2 = imp_line.split('         ')
            #Print and check ME filling
            MG_ggnoHZZ_ME_value = float(ME_list2[1])
            #Remove the files created
            os.remove(results)
            #going back to home directory to start another calculation
            os.chdir(TarArea)

            ########################
            #   NLO qqZZ ME-inc    #
            ########################

            #Going to the right directory
            if(mytree.event_type<2):
                os.chdir(ME_qqZZ_4l)
                subprocess.call(["./check",ps_input,results])
            if(mytree.event_type>=2):
                os.chdir(ME_qqZZ_2l2l)
                subprocess.call(["./check",ps_input,results])
            #Read the ME from the results file
            file = open(results, "r")
            outputLines = file.readlines()
            #Read the ME from the results file
            imp_line_1 = outputLines[7]
            imp_line_2 = outputLines[22]
            ME_born = imp_line_1.split('         ')
            ME_finite = imp_line_2.split('    ')
            #Print and check ME filling
            MG_NLO_qqZZ_ME_value =( float(ME_born[1])+ float(ME_finite[1]))

            #Remove the files created
            os.remove(results)
            #going back to home directory to start another calculation
            os.chdir(TarArea)
            
            #filling the 0 jet events with just 0 jet ME and 99 for 1 jet MEs
            MG_ggHZZ_ME[0] = MG_ggHZZ_ME_value
            MG_ggnoHZZ_ME[0] = MG_ggnoHZZ_ME_value
            MG_NLO_qqZZ_ME[0] = MG_NLO_qqZZ_ME_value


            newTree.Fill()
            
            os.chdir(TarArea)                   
        #Output the different ME for check
        print("#####################################################################################################################################################")
        print(newTree.MG_ggHZZ_ME)
        print(newTree.MG_ggnoHZZ_ME)
        print(newTree.MG_NLO_qqZZ_ME)

    
        os.chdir(TarArea)
        os.chdir(scratchArea)

    f.Write()
    f.Close()

    file_1 = "tmp_"+process+"_"+str(i)
    subprocess.call(["rm","-rf",file_1])
    #os.remove(file_1)

    minitree = "./"+output_tree
    #subprocess.call("pwd")
    #subprocess.call("ls")

    #copy the tree back to gpfs
    subprocess.call(["scp", minitree ,"abc-at13:/gpfs3/umass/HZZ/MG_0j_Ana_Offshell/post_process_minitrees/"])
    subprocess.call(["rm","-rf",output_tree])

    print (datetime.now()-start)
    
if __name__ == "__main__":
    import sys
    Calculate_ME(int(sys.argv[1]),int(sys.argv[2]),str(sys.argv[3]),str(sys.argv[4]))


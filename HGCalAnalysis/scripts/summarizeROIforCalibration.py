import ROOT
import array
import sys

from RecoNtuples.HGCalAnalysis.ROISummary import ROISummary


def main(url_in=sys.argv[1],url_out=sys.argv[2]):
    
    """makes a simpler summary of the ROI for fast calibration"""
    
    # open the output file
    fout = ROOT.TFile(url_out, "RECREATE")
    fout.cd()
    varnames  = ['genen1','geneta1','genphi1']
    varnames += ['genen2','geneta2','genphi2']
    for i in xrange(1,4):
        varnames+=['en1_%d'%i,'noise1_%d'%i,'en2_%d'%i,'noise2_%d'%i]
    output_tuple = ROOT.TNtuple("data","data",":".join(varnames))

    #loop over events
    fIn=ROOT.TFile.Open(url_in)
    t=fIn.Get('ana/data')
    for i in xrange(0,t.GetEntriesFast()):

        t.GetEntry(i)

        #define the ROIs
        roiList={}
        for ir in xrange(0,t.ROIs.size()):
            if t.ROIs[ir].id()<0 : continue
            roiList[ir]=ROISummary(t.ROIs[ir].p4())

        for h in t.Hits:
            roiIdx  = h.associatedROI()
            rid     = t.ROIs[roiIdx].id()
            roiKey  = roiIdx if rid>0 else abs(rid)-1
            en      = h.en()
            isNoise = True if rid<0 else False
            regIdx  = h.signalRegion()
            roiList[roiKey].addHit(en=en,isNoise=isNoise,regIdx=regIdx)
            
        #fill this event in the ntuple
        varvals=[]

        #mc truth
        for ir in xrange(0,2):
            genP4=ROOT.TLorentzVector(0,0,0,0)
            if ir in roiList:
                genP4=roiList[ir].genP4
            varvals += [genP4.E(),genP4.Eta(),genP4.Phi()]
            
        #results for the 3 different regions
        for ireg in xrange(1,4):
            for ir in xrange(0,2):
                recP4=ROOT.TLorentzVector(0,0,0,0)
                noise=0.
                if ir in roiList:
                    recP4    = roiList[ir].getReconstructedP4(ireg)
                    noise    = roiList[ir].getNoiseInRing(ireg)/5. #averaged over 5 different regions
                varvals += [recP4.E(),noise]        

        #all done
        output_tuple.Fill(array.array("f",varvals))

    #close file
    fout.cd()
    output_tuple.Write()
    fout.Close()


if __name__ == "__main__":
    main()
            

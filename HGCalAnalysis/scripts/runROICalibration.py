import ROOT
import numpy as np
from array import array
import sys
from collections import OrderedDict
from RecoNtuples.HGCalAnalysis.RootTools import buildMedianProfile

MINGENEN=100
MINGENETA=1.6
MAXGENETA=2.8
A=[ROOT.TMath.Pi()*1.3**2, ROOT.TMath.Pi()*2.6**2, ROOT.TMath.Pi()*5.3**2]

def getEnergiesForCalibration(data,ireg):

    """returns an array of reconstructed energies, eta and average noise"""

    x=[]
    for i in xrange(0,data.GetEntriesFast()):

        data.GetEntry(i)

        for ia in xrange(1,3):
            genen=getattr(data,'genen%d'%ia)
            geneta=abs(getattr(data,'geneta%d'%ia))
            if genen<MINGENEN: continue
            if geneta<MINGENETA: continue
            if geneta>MAXGENETA: continue
            recen=getattr(data,'en%d_%d'%(ia,ireg))
            #avgnoise=getattr(data,'noise%d_%d'%(ia,ireg))
            avgnoise=getattr(data,'noise%d_3'%ia)*A[ireg-1]/A[2]
            x.append( [genen,geneta,recen,avgnoise] )
  

    return np.array(x)

def calibrateSpectrum(h,title,proc,func='pol1'):

    """calibrates a DeltaE/E versus x spectrum based on the profile"""
    
    c=ROOT.TCanvas('c','c',500,500)
    c.SetTopMargin(0.05)
    c.SetBottomMargin(0.1)
    c.SetLeftMargin(0.12)
    c.SetRightMargin(0.1)
    #prof=h.ProfileX()
    prof=buildMedianProfile(h)
    h.Draw('colz')
    prof.Draw('e2p')
    prof.Fit(func)
    calibGr=prof.GetListOfFunctions().At(0).Clone('%s%s_calib'%(h.GetName(),title))        
    tex=ROOT.TLatex()
    tex.SetTextFont(42)
    tex.SetTextSize(0.04)
    tex.SetNDC()
    tex.DrawLatex(0.12,0.96,'#bf{CMS} #it{simulation preliminary}')
    tex.DrawLatex(0.15,0.88,title)
    tex.SetTextAlign(31)
    tex.DrawLatex(0.97,0.96,proc)        
    c.SaveAs('%s%s.png'%(h.GetName(),title))
    prof.Delete()        
    c.Delete()

    return calibGr

def doL0L1Calibration(url,calib,nq=6):

    """performs the L0 (uniform eta response) and L1 (absolute scale) calibrations"""

    #derive L0 and L1 calibrations for 3 different signal regions
    fIn=ROOT.TFile.Open(url)
    data=fIn.Get('data')
    for ireg in [1,2,3]:
        x=getEnergiesForCalibration(data,ireg)
        xq=np.percentile(x, [i*100./nq for i in xrange(0,nq+1)], axis=0)

        #relative calibration versus eta        
        resolVsEta = ROOT.TH2F('resolvseta',';|#eta|;#DeltaE/E', nq, array('d',xq[:,1]), 100,-1,1)
        for i in xrange(0,len(x)):
            genEn, genEta, recEn,_ = x[i]
            deltaE=recEn/genEn-1.
            resolVsEta.Fill(genEta,deltaE)
        calib['L0'][ireg]=calibrateSpectrum(resolVsEta,'SR%d'%ireg,'H#rightarrow#gamma#gamma (PU=0)','pol2')
    
        #relative calibration versus energy
        resolVsEn = ROOT.TH2F('resolvsen',';Reconstructed energy [GeV];#DeltaE/E', nq, array('d',xq[:,0]), 100,-1,1)
        for i in xrange(0,len(x)):
            genEn,genEta,recEn,_=x[i]
            recEn=recEn/(calib['L0'][ireg].Eval(genEta)+1.0)
            deltaE=recEn/genEn-1.
            resolVsEn.Fill(recEn,deltaE)
            resolVsEta.Fill(genEta,deltaE)
        calib['L1'][ireg]=calibrateSpectrum(resolVsEn,'SR%d'%ireg,'H#rightarrow#gamma#gamma (PU=0)','pol1')

    fIn.Close()

    return calib

def doPUCalibration(url,calib,nq=10):

    """parametrizes the absolute shift in energy as function of the average noise in the SR cone"""

    fIn=ROOT.TFile.Open(url)
    data=fIn.Get('data')

    for ireg in [1,2,3]:
        x=getEnergiesForCalibration(data,ireg)
        xq=np.percentile(x, [i*100./nq for i in xrange(0,nq+1)], axis=0)

        #relative calibration versus eta
        resolVsNoiseFrac = ROOT.TH2F('resolvsnoisefrac',';<Noise> [GeV];#DeltaE [GeV]', nq,array('d',xq[:,3]),50,-100,100)

        for i in xrange(0,len(x)):
            genEn,genEta,recEn,noiseEst=x[i]
            recEn=recEn/(calib['L0'][ireg].Eval(genEta)+1.0)
            recEn=recEn/(calib['L1'][ireg].Eval(recEn)+1.0)
            deltaE=recEn-genEn 
            resolVsNoiseFrac.Fill(noiseEst,deltaE)
    
        calib['L2'][ireg]=calibrateSpectrum(resolVsNoiseFrac,'SR%d'%ireg,'H#rightarrow#gamma#gamma (PU=140)','pol2')

    return calib


def applyCalibrationTo(url,calib,title,silent=False):
    
    """applies the calibration to the photons and shows the energy and H->gg mass resolutions """
    
    pfix=''.join(calib.keys())


    fIn=ROOT.TFile.Open(url)
    data=fIn.Get('data')

    histos={}
    for ireg in xrange(1,4):
        histos['dm%d'%ireg]  = ROOT.TH1F('dm%d'%ireg, ';#Delta m_{#gamma#gamma}/m_{#gamma#gamma};PDF',50,-0.1,0.1)
        histos['den%d'%ireg] = ROOT.TH1F('den%d'%ireg,';#Delta E/E;PDF',50,-0.1,0.1)
    for h in histos:
        histos[h].Sumw2()
        histos[h].SetLineColor(1)
        histos[h].SetMarkerColor(1)
        histos[h].SetMarkerStyle(20)
        histos[h].SetDirectory(0)

    for i in xrange(0,data.GetEntriesFast()):
        
        data.GetEntry(i)

        #generator level photons
        genphotons=[]
        for ia in xrange(1,3):
            genen  = getattr(data,'genen%d'%ia)
            geneta = getattr(data,'geneta%d'%ia)
            genphi = getattr(data,'genphi%d'%ia)
            genphotons.append(ROOT.TLorentzVector(0,0,0,0))
            genphotons[-1].SetPtEtaPhiM(genen/ROOT.TMath.CosH(geneta),geneta,genphi,0.)
        genh = genphotons[0]+genphotons[1]

        #H->gg fiducial cuts
        if genphotons[0].Pt()<20 or  genphotons[1].Pt()<20 : continue
        if genphotons[0].Pt()<40 and genphotons[1].Pt()<40 : continue
        if abs(genphotons[0].Eta())<1.5 or abs(genphotons[1].Eta())<1.5 : continue
        if abs(genphotons[0].Eta())>2.8 or abs(genphotons[1].Eta())>2.8 : continue

        #reconstructed photons in different regions
        for ireg in xrange(1,4):

            photons=[]
            for ia in xrange(1,3):
                genen    = getattr(data,'genen%d'%ia)
                if genen<MINGENEN : continue
                geneta   = getattr(data,'geneta%d'%ia)
                genphi   = getattr(data,'genphi%d'%ia)
                recen    = getattr(data,'en%d_%d'%(ia,ireg))
                #avgnoise = getattr(data,'noise%d_%d'%(ia,ireg))
                avgnoise=getattr(data,'noise%d_3'%ia)*A[ireg-1]/A[2]

                if 'L0' in calib:
                    recen=recen/(calib['L0'][ireg].Eval(abs(geneta))+1.0)
                    if 'L1' in calib:
                        recen=recen/(calib['L1'][ireg].Eval(recen)+1.0)
                        if 'L2' in calib and ireg in calib['L2']:
                            recen=recen-calib['L2'][ireg].Eval(avgnoise)

                deltaE = recen/genen-1.
                histos['den%d'%ireg].Fill(deltaE)
                photons.append(ROOT.TLorentzVector(0,0,0,0))
                photons[-1].SetPtEtaPhiM(recen/ROOT.TMath.CosH(geneta),geneta,genphi,0.)
            
            if len(photons)!=2: continue
            h = photons[0]+photons[1]
            deltaM=h.M()/genh.M()-1
            histos['dm%d'%ireg].Fill(deltaM)

    fIn.Close()

    if silent: return histos

    c=ROOT.TCanvas('c','c',500,500)
    c.SetTopMargin(0.05)
    c.SetBottomMargin(0.1)
    c.SetLeftMargin(0.12)
    c.SetRightMargin(0.03)
    for ireg in xrange(1,4):
        for k in ['dm','den']:
            h=histos['%s%d'%(k,ireg)]
            h.Scale(1./h.Integral())
            h.Draw()
            h.GetYaxis().SetTitleOffset(0.9)
            h.GetYaxis().SetRangeUser(0,h.GetMaximum()*1.2)
            h.Fit('gaus','M+')
            gaus=h.GetListOfFunctions().At(0)
            tex=ROOT.TLatex()
            tex.SetTextFont(42)
            tex.SetTextSize(0.04)
            tex.SetNDC()
            tex.DrawLatex(0.12,0.96,'#bf{CMS} #it{simulation preliminary}')
            tex.DrawLatex(0.15,0.88,'SR%d (%s-calibrated)'%(ireg,pfix))
            tex.DrawLatex(0.15,0.84,'#mu=%3.3f#pm%3.3f'%(gaus.GetParameter(1),gaus.GetParError(1)))
            tex.DrawLatex(0.15,0.80,'#sigma=%3.3f#pm%3.3f'%(gaus.GetParameter(2),gaus.GetParError(2)))
            tex.SetTextAlign(31)
            tex.DrawLatex(0.97,0.96,title)
            c.SaveAs('%s%s.png'%(pfix,h.GetName()))

    #save in a local file
    fOut=ROOT.TFile.Open('calib%s.root'%pfix,'RECREATE')
    for h in histos: histos[h].Clone().Write()
    fOut.Close()

    return histos

def format(pList,name,ci,fill=-1,fillci=None):
    for p in pList:
        pList[p].SetName('%s%s'%(p,name))
        pList[p].SetLineWidth(2)
        pList[p].SetLineColor(ci)
        pList[p].Scale(1./pList[p].Integral())
        if fill<0: continue
        pList[p].SetFillStyle(fill)
        if not fillci: fillci=ci
        pList[p].SetFillColor(fillci)

def main():

    ROOT.gStyle.SetOptStat(0)
    ROOT.gStyle.SetOptTitle(0)
    ROOT.gROOT.SetBatch(True)
    ROOT.gStyle.SetPalette(ROOT.kTemperatureMap)

    toCompare=OrderedDict()
    
    nopuF=sys.argv[1]
    calib=OrderedDict()
    calib['L0']={}
    calib['L1']={}
    doL0L1Calibration(url=nopuF,calib=calib)
    toCompare['PU=0']=applyCalibrationTo(url=nopuF,calib=calib,title='H#rightarrow#gamma#gamma (PU=0)')
    format(toCompare['PU=0'],'nopu',1,1001,ROOT.kGray)

    puF=sys.argv[2]
    puTag=sys.argv[3]
    calib['L2']={}
    doPUCalibration(url=puF,calib=calib)
    toCompare['PU=%s'%puTag]=applyCalibrationTo(url=puF,calib=calib,title='H#rightarrow#gamma#gamma (PU=%s)'%puTag)
    format(toCompare['PU=%s'%puTag],puTag,ROOT.kRed+2)

    #save final calibration
    import pickle
    with open('calib_pu%s.pck'%puTag,'w') as cachefile:
        pickle.dump(calib,cachefile, pickle.HIGHEST_PROTOCOL)
   
    if len(sys.argv)>4:
        iextra=0
        extraColors=[ROOT.kAzure+7, ROOT.kBlue-7,ROOT.kMagenta, ROOT.kMagenta+2, ROOT.kMagenta-9]
        for x in sys.argv[4].split(','):
            tag,url=x.split(':')
            toCompare[tag]=applyCalibrationTo(url=url,calib=calib,title='H#rightarrow#gamma#gamma (%s)'%tag,silent=True)
            format(toCompare[tag],'extra%d'%iextra,extraColors[iextra])
            iextra+=1

    showFinalComparison(toCompare,title='H#rightarrow#gamma#gamma (PU=%s)'%puTag)


def showFinalComparison(compList,title):
    c=ROOT.TCanvas('c','c',500,500)
    c.SetTopMargin(0.05)
    c.SetBottomMargin(0.1)
    c.SetLeftMargin(0.12)
    c.SetRightMargin(0.03)
    plots=compList[compList.keys()[0]].keys()
    for p in plots:
        c.Clear()

        leg=ROOT.TLegend(0.65,0.9,0.95,0.7)
        leg.SetBorderSize(0)
        leg.SetFillStyle(0)
        leg.SetTextFont(42)

        drawOpt='hist'
        for tag in compList:
            compList[tag][p].Rebin()
            compList[tag][p].Draw(drawOpt)
            drawOpt='histsame'
            compList[tag][p].GetYaxis().SetTitleOffset(0.9)
            compList[tag][p].GetYaxis().SetRangeUser(0,compList[tag][p].GetMaximum()*1.2)
            leg.AddEntry(compList[tag][p],tag,'lf')

        leg.Draw()

        tex=ROOT.TLatex()
        tex.SetTextFont(42)
        tex.SetTextSize(0.04)
        tex.SetNDC()
        tex.DrawLatex(0.12,0.96,'#bf{CMS} #it{simulation preliminary}')
        tex.SetTextAlign(31)
        tex.DrawLatex(0.97,0.96,title)
        c.SaveAs('comp_%s.png'%(compList[tag][p].GetName()))


if __name__ == "__main__":
    main()
              

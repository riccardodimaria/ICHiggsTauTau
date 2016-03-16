#!/usr/bin/env python
import ROOT
import CombineHarvester.CombineTools.plotting as plot
from UserCode.ICHiggsTauTau.uncertainties import ufloat


def BuildVarString(var):
    pos = var.rfind('(')
    newvar = var[:pos] + '>>htemp' + var[pos:]
    return newvar


def BuildCutString(sel, weight):
    return '(%s)&&(%s)' % (sel, weight)


def Rate(hist):
    err = ROOT.Double()
    integral = hist.IntegralAndError(0, hist.GetNbinsX()+1, err)
    return ufloat(integral, err)


def Shape(tree, var, sel='1', weight='1'):
    ROOT.TH1.SetDefaultSumw2(True)
    # print BuildVarString(var)
    # print BuildCutString(sel, weight)
    ROOT.TH1.AddDirectory(True)
    tree.Draw(BuildVarString(var),  BuildCutString(sel, weight), 'goff')
    ROOT.TH1.AddDirectory(False)
    res = ROOT.gDirectory.Get('htemp').Clone()
    ROOT.gDirectory.Delete('htemp;*')
    return res


def AsTex(val):
    return '%.2f#pm%.2f' % (val.nominal_value, val.std_dev)


## Boilerplate
ROOT.PyConfig.IgnoreCommandLineOptions = True
ROOT.gROOT.SetBatch(ROOT.kTRUE)
plot.ModTDRStyle()
# ROOT.gStyle.SetNdivisions(510, 'XYZ') # probably looks better

for sel in [
        ('iso_m_iso_t', 'iso_1<0.1 && iso_2<1.5'),
        ('anti_m_iso_t', 'iso_1>0.2 && iso_2<1.5'),
        ('anti_m_iso_t_btag1', 'iso_1>0.2 && iso_2<1.5 && nbtags>=1'),
        ('anti_m_anti_t', 'iso_1>0.2 && iso_2>1.5'),
        ('iso_m_anti_t', 'iso_1<0.1 && iso_2>1.5')]:

    for draw in [
            ('TMath::Abs(flavour_1)(25,-0.5,24.5)', ROOT.kBlue, 'Muon partonFlavor'),
            ('TMath::Abs(flavour_2)(25,-0.5,24.5)', ROOT.kRed, 'Tau partonFlavor')]:

        draw_var = draw[0]

        name = draw_var[:draw_var.rfind('(')]
        name = name.replace('::', '.')
        name = name.replace('(', '_')
        name = name.replace(')', '_')

        canv = ROOT.TCanvas(sel[0]+'__'+name, sel[0]+'__'+name)
        pads = plot.TwoPadSplit(0.30, 0.01, 0.01)
        # pads = plot.OnePad()

        file = ROOT.TFile('QCD.root')
        tree = file.Get('qcdstudy')

        # How to make a discrete map?

        hist_ss = Shape(tree, draw_var, '(q_1 == q_2) && (%s)' % sel[1])
        rate_ss = Rate(hist_ss)
        hist_os = Shape(tree, draw_var, '(q_1 != q_2) && (%s)' % sel[1])
        rate_os = Rate(hist_os)

        plot.Set(hist_ss, LineWidth=2, LineColor=14, MarkerColor=14)
        plot.Set(hist_os, LineWidth=2, LineColor=draw[1], MarkerColor=draw[1])

        print 'rate_os: ' + str(rate_os)
        print 'rate_ss: ' + str(rate_ss)

        os_ss_ratio = rate_os / rate_ss
        print 'ratio: ' + str(os_ss_ratio)
        hist_os.Draw('HISTE')
        hist_ss.Draw('HISTSAME')
        hist_ss.Draw('ESAME')
        # plot.DrawCMSLogo(pads[0], 'CMS', args.cms_sub, 11, 0.045, 0.035, 1.2, '', 0.8)
        # plot.DrawTitle(pads[0], args.title_right, 3)
        # plot.DrawTitle(pads[0], args.title_left, 1)
        plot.FixTopRange(pads[0], plot.GetPadYMax(pads[0]), 0.3)

        latex = ROOT.TLatex()
        latex.SetNDC()
        latex.SetTextSize(0.04)
        latex.DrawLatex(0.2, 0.85, "OS/SS = " + AsTex(os_ss_ratio))

        axis = plot.GetAxisHist(pads[0])
        plot.Set(axis.GetXaxis(), Title=draw[2])
        plot.Set(axis.GetYaxis(), Title='Events')

        legend = ROOT.TLegend(0.65, 0.75, 0.87, 0.90, '', 'NBNDC')

        legend.AddEntry(hist_os, 'OS events', 'LF')
        legend.AddEntry(hist_ss, 'SS events', 'LF')
        legend.Draw()

        pads[1].cd()
        pads[1].SetGrid(0, 1)

        ratio = hist_os.Clone()
        ratio.Divide(hist_ss)
        ratio.Draw()


        plot.SetupTwoPadSplitAsRatio(pads, plot.GetAxisHist(pads[0]), plot.GetAxisHist(pads[1]), 'OS/SS', True, 0.05, 1.95)

        # plot.DrawTitle(pads[0], "OS/SS = " + AsTex(os_ss_ratio), 1)
        plot.DrawTitle(pads[0], 'QCD_Pt_20_MuEnrichedPt_15, 8 TeV', 3)

        canv.Print('.pdf')
        canv.Print('.png')
# maketable.TablefromJson(args.table_vals, args.file, "TablefromJson.txt")

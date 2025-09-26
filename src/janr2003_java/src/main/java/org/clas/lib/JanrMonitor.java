package org.clas.lib;

import java.awt.BorderLayout;
import java.awt.Dimension;
import java.awt.FlowLayout;
import java.awt.FontMetrics;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.Insets;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.ItemEvent;
import java.awt.event.ItemListener;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.Set;

import javax.swing.BorderFactory;
import javax.swing.ButtonGroup;
import javax.swing.DefaultComboBoxModel;
import javax.swing.JButton;
import javax.swing.JCheckBox;
import javax.swing.JComboBox;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JRadioButton;
import javax.swing.JSlider;
import javax.swing.JSplitPane;
import javax.swing.SwingConstants;
import javax.swing.SwingUtilities;
import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;

import org.jlab.groot.base.GStyle;

public class JanrMonitor implements ActionListener {
	
	Janr03 j = new Janr03();

    private JPanel                     engineView = new JPanel();
    private EmbeddedCanvasTabbed             data = new EmbeddedCanvasTabbed("DATA");
    private EmbeddedCanvasTabbed             fitc = new EmbeddedCanvasTabbed("FITC");

    private String                   detectorName = null;
    private JPanel                  detectorPanel = new JPanel();
    
    private JSplitPane                      vPane = null;
    
    private EmbeddedCanvasTabbed[] detectorCanvas = new EmbeddedCanvasTabbed[2];
    private ArrayList<String>[]  detectorTabNames = (ArrayList<String>[])new ArrayList[2];
    public  DataGroupManager[]                dgm = new DataGroupManager[2];
      
    public Boolean                    dgmActive = false;
    public int                        runNumber = 0;

    private int                         bin1234 = 0;
    private int                          binWCF = 0;
    public  int                             nt  = 0;
    private int                             bin = 0;
    public int                              cof = 0;
    public Set<Integer>                 cofList = new HashSet<>();
    		
    private JPanel                  actionPanel = null;
    private JPanel               controlsPanel0 = null;
    private JPanel               controlsPanel1 = null;
 
    private JComboBox                 configCMB = new JComboBox();
	private GridBagConstraints              gbc = new GridBagConstraints();

    public Boolean                  histosExist = false;   
    public String                          root = " ";
    public JSlider                    binslider = null;
    public JSlider                    cofslider = new JSlider(JSlider.HORIZONTAL,10, 100, 50); 
    private ButtonGroup bG0, bG1;
    private JRadioButton bW,bC,bF,b1,b2,b3,b4;
    private Boolean useWCF=false, use1234=false, useBSP=false;

    private JLabel  cof_resonance = new JLabel();
    private JLabel cof_parameters = new JLabel();
    private JLabel      cof_label = new JLabel("" + String.format("%d", 0)); 
    
    public JanrMonitor(String name) {
		initGStyle(14);
	    detectorName   = name;	    
        root           = detectorName+".JanrMonitor.";	
    	for (int i=0; i<2; i++) detectorTabNames[i] = new ArrayList<String>();
    	for (int i=0; i<2; i++) dgm[i] = new DataGroupManager();
    	for (int i=0; i<2; i++) detectorCanvas[i] = new EmbeddedCanvasTabbed();
	}
    
    public void init() {
    	System.out.println(root+"init");
        initPanel();    	
    }
	
    public void initGStyle(int fontsize) {
        GStyle.getAxisAttributesX().setTitleFontSize(fontsize);
        GStyle.getAxisAttributesX().setLabelFontSize(fontsize);
        GStyle.getAxisAttributesY().setTitleFontSize(fontsize);
        GStyle.getAxisAttributesY().setLabelFontSize(fontsize);
        GStyle.getAxisAttributesZ().setLabelFontSize(fontsize); 
        GStyle.getAxisAttributesX().setAxisGrid(false);
        GStyle.getAxisAttributesY().setAxisGrid(false);
        GStyle.getAxisAttributesX().setLabelFontName("Avenir");
        GStyle.getAxisAttributesY().setLabelFontName("Avenir");
        GStyle.getAxisAttributesZ().setLabelFontName("Avenir");
        GStyle.getAxisAttributesX().setTitleFontName("Avenir");
        GStyle.getAxisAttributesY().setTitleFontName("Avenir");
        GStyle.getAxisAttributesZ().setTitleFontName("Avenir");
        GStyle.getAxisAttributesZ().setAxisAutoScale(true);    	
        GStyle.getGraphErrorsAttributes().setMarkerStyle(1);
        GStyle.getGraphErrorsAttributes().setMarkerColor(2);
        GStyle.getGraphErrorsAttributes().setMarkerSize(2);
        GStyle.getGraphErrorsAttributes().setLineColor(2);
        GStyle.getGraphErrorsAttributes().setLineWidth(1);
        GStyle.getGraphErrorsAttributes().setFillStyle(1);   
    }	
    
    public JPanel getPanel(EmbeddedCanvasTabbed p1, EmbeddedCanvasTabbed p2) {        
        engineView.setLayout(new BorderLayout());
        engineView.add(getCanvasPane(p1,p2),BorderLayout.CENTER);
        return engineView;       
     }  
    
     public JSplitPane getCanvasPane(EmbeddedCanvasTabbed p1, EmbeddedCanvasTabbed p2) {
        p1.setPreferredSize(new Dimension(2300,1000)); p2.setPreferredSize(new Dimension(2300,100));
        vPane = new JSplitPane(JSplitPane.VERTICAL_SPLIT,p1,p2);    	  
        vPane.setDividerLocation(0.8);
        return vPane;
     }
    
    public void initPanel() {
        detectorPanel.setLayout(new BorderLayout());
        actionPanel    = new JPanel(new FlowLayout());
        controlsPanel0 = new JPanel(new GridBagLayout());        
        controlsPanel1 = new JPanel();
               
        gbc.fill = GridBagConstraints.HORIZONTAL; 
        gbc.weightx = 0.5; gbc.gridx=0 ; gbc.gridy=0 ;        
        controlsPanel0.add(controlsPanel1,gbc);
        
        controlsPanel1.setBorder(BorderFactory.createTitledBorder("Bins"));		       
        controlsPanel1.add(packWCFPanel());
        controlsPanel1.add(packBinSliderPanel());		       
        controlsPanel1.add(packComboBox());
        controlsPanel1.add(packCofResonance());
        controlsPanel1.add(packCofSliderPanel());
        controlsPanel1.add(packCofSliderLabelPanel());
        controlsPanel1.add(packCofParameters());
        controlsPanel1.add(packCofParametersLabel());
        controlsPanel1.add(packResetButton());

//        getDetectorPanel().add(getPanel(data,fitc),BorderLayout.CENTER);           
        detectorPanel.add(getPanel(getDetectorCanvas(0),getDetectorCanvas(1)),BorderLayout.CENTER);  //DATA, FITC canvases         
//        getDetectorPanel().add(getDetectorCanvas,BorderLayout.CENTER);           
        detectorPanel.add(controlsPanel0,     BorderLayout.SOUTH); 
        
        vPane.setResizeWeight(0.75);

    }
    
    public JPanel packWCFPanel() {
    	if (useWCF || use1234) actionPanel.add(getWCFButtonPanel());
    	return actionPanel;
    }
    
    public JPanel packBinSliderPanel() {
        actionPanel.add(getBinSliderPanel());	
        return actionPanel;
    }
    
    public JPanel packComboBox() {
    	actionPanel.add(getComboBox());
    	return actionPanel;
    }

    public JPanel packCofResonance() {
    	actionPanel.add(getCofResonance());
    	return actionPanel;
    }
          
    public JPanel packCofSliderPanel() {
    	actionPanel.add(getCofSliderPanel());
    	return actionPanel;
    }
    
    public JPanel packCofSliderLabelPanel() {
     	actionPanel.add(getCofSliderLabelPanel());
    	return actionPanel;
    } 
    
    public JPanel packCofParameters() {
    	actionPanel.add(getCofParameters());
    	return actionPanel;
    } 
    
    public JPanel packCofParametersLabel() {
    	actionPanel.add(getCofParmametersLabel());
    	return actionPanel;
    } 
    
    public JPanel packResetButton() {
    	actionPanel.add(getResetButton());
    	return actionPanel;
    } 

    public void use1234Buttons(boolean flag) {
    	use1234 = flag;
    }    
    
    public void useWCFButtons(boolean flag) {
    	useWCF = flag;
    }    
    
    public void useBinSliderPane(boolean flag) {
    	useBSP = flag;
    }
    
    public JPanel getWCFButtonPanel() {
        JPanel buttonPane = new JPanel();         
        if(useWCF) {
        	bW = new JRadioButton("W"); buttonPane.add(bW); bW.setActionCommand("0"); bW.addActionListener(this);
        	bC = new JRadioButton("C"); buttonPane.add(bC); bC.setActionCommand("1"); bC.addActionListener(this); 
        	bF = new JRadioButton("F"); buttonPane.add(bF); bF.setActionCommand("2"); bF.addActionListener(this);
        	bG0 = new ButtonGroup(); bG0.add(bW); bG0.add(bC); bG0.add(bF); 
        	bW.setSelected(true); bW.doClick();
        }
        if(use1234) {
            b1 = new JRadioButton("1"); buttonPane.add(b1); b1.setActionCommand("0"); b1.addActionListener(this); 
            b2 = new JRadioButton("2"); buttonPane.add(b2); b2.setActionCommand("1"); b2.addActionListener(this); 
            b3 = new JRadioButton("3"); buttonPane.add(b3); b3.setActionCommand("2"); b3.addActionListener(this); 
            b4 = new JRadioButton("4"); buttonPane.add(b4); b4.setActionCommand("3"); b4.addActionListener(this);
            bG1 = new ButtonGroup(); bG1.add(b1); bG1.add(b2); bG1.add(b3); bG1.add(b4);        	
        	b1.setSelected(true); b1.doClick();
        }
        return buttonPane;
    }
    
    public JPanel getBinSliderPanel() {
        JPanel sliderPane = new JPanel();
        JLabel       label = new JLabel("" + String.format("%d", 0));       
        binslider = new JSlider(JSlider.HORIZONTAL, 0, 30, 10); 
        binslider.setPreferredSize(new Dimension(100,10));
        sliderPane.add(new JLabel("WCF:",JLabel.CENTER));
        sliderPane.add(binslider);
        sliderPane.add(label);  
        binslider.addChangeListener(new ChangeListener() {
            public void stateChanged(ChangeEvent e) {            	
                JSlider slider1 = (JSlider) e.getSource(); 
                bin = slider1.getValue();           
                label.setText(String.valueOf(""+String.format("%d", bin)));
                plotHistos(getRunNumber(), nt, bin);
            }
        });         
        return sliderPane;
    }
   
    public JPanel getComboBox() {
		JPanel sliderParm = new JPanel();
		sliderParm.setLayout(new FlowLayout());      
		sliderParm.add(new JLabel("COEF:"));  				
		DefaultComboBoxModel model = (DefaultComboBoxModel) configCMB.getModel();
		for (int i=0; i<j.pname.length; i++) model.addElement(j.pname[i]);
		configCMB.setModel(model);
		configCMB.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				cof = model.getIndexOf(configCMB.getSelectedItem()); 
				cofList.add(cof); 
				cof_resonance.setText(j.res_name[cof]);
				cofslider.setValue(((int)j.start_value[cof]+5)*10); //new
				Res2Amp(j.res_name[cof]);
                plotHistos(getRunNumber(), nt, bin);
			}
		});
		cofList.add(0);
        configCMB.setSelectedIndex(0);
		sliderParm.add(configCMB); 
		
		return sliderParm;
    }

    public JPanel getCofSliderPanel() {
        JPanel sliderPanel = new JPanel(new GridBagLayout());                     
        cofslider.addChangeListener(new ChangeListener() {
            public void stateChanged(ChangeEvent e) {                                 			
                j.xnew[cof] = cofslider.getValue()*0.1f-5f;  
                cof_label.setText(String.valueOf(""+String.format("%.2f", j.xnew[cof])));
                hjanr_loadpar(2);
                Res2Amp(j.res_name[cof]);
                plotHistos(getRunNumber(), nt, bin);
            }
        });    
        cofslider.setValue(((int)j.start_value[0]+5)*10); //default for startup
        gbc.gridx=0; gbc.gridy=0;
        sliderPanel.add(cofslider,gbc);
        return sliderPanel;
    }
    
    public JPanel getCofSliderLabelPanel() {
    	JPanel labelPanel = new JPanel(new GridBagLayout());
    	
        // Get the widest possible text (to avoid slider jitter when sign changes)
        String widestText = "-4.00";
        FontMetrics fm = cof_label.getFontMetrics(cof_label.getFont());
        int widestWidth = fm.stringWidth(widestText);
        int widestHeight = fm.getHeight();
        
        // Set the preferred size to accommodate the widest string
        Dimension preferredSize = new Dimension(widestWidth, widestHeight);
        cof_label.setPreferredSize(preferredSize);
        
        gbc.insets = new Insets(10,10,10,10);
        gbc.gridx=1; gbc.gridy=0;
    	labelPanel.add(cof_label,gbc);
    	return labelPanel;
    }
    
    public JPanel getCofParmametersLabel() {
    	JPanel labelPanel = new JPanel(new GridBagLayout());
    	
        // Get the widest possible text (to avoid slider jitter when sign changes)
        String widestText = "xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx";
        FontMetrics fm = cof_parameters.getFontMetrics(cof_parameters.getFont());
        int widestWidth = fm.stringWidth(widestText);
        int widestHeight = 2*fm.getHeight();
        
        // Set the preferred size to accommodate the widest string
        Dimension preferredSize = new Dimension(widestWidth, widestHeight);
        cof_parameters.setPreferredSize(preferredSize);
        
        gbc.insets = new Insets(1,1,1,1);
        gbc.gridx=1; gbc.gridy=0;
    	labelPanel.add(cof_parameters,gbc);
    	return labelPanel;
    }
        
    public JLabel getCofParameters() {
        return cof_parameters;
    }
    
    public JLabel getCofResonance() {
        return cof_resonance;
    }

    public void Res2Amp(String res) {
    	double am=0, ae=0, as=0, a32=0, a12=0, s12=0, rem=0, rsm=0;
    	hjanr_loadpar(2); j.jip.janrIniPoint(0.4);
    	switch (res) {
    	case "P33(1232)": am=j.am3[0];    ae=j.ae3[0];    as=j.as3[0]; rem=ae/am; rsm=as/am; break;
    	case "P11(1440)": am=j.am1[0][0]; ae=0;           as=j.as1[0][0]; a32=0;                a12=j.aa1[21]*j.cmp1; s12=j.sa1[21]*j.csp1; break;
    	case "S11(1535)": am=0;           ae=j.ae1[0][1]; as=j.as1[0][1]; a32=0;                a12=j.aa1[22]*j.cep2; s12=j.sa1[22]*j.csp2; break;
    	case "D13(1520)": am=j.am1[0][2]; ae=j.ae1[0][2]; as=j.as1[0][2]; a32=j.aa3[23]*j.cep3; a12=j.aa1[23]*j.cmp3; s12=j.sa1[23]*j.csp3; break;
    	case "S11(1650)": am=0;           ae=j.ae1[0][4]; as=j.as1[0][4]; a32=0;                a12=j.aa1[25]*j.cep4; s12=j.sa1[25]*j.csp4; break;
    	case "S31(1620)": am=0;           ae=j.ae3[4];    as=j.as3[4];    a32=0;                a12=j.aa1[ 5]*j.cep5; s12=j.sa1[ 5]*j.csp5; break;
    	case "F15(1680)": am=j.am1[0][5]; ae=j.ae1[0][5]; as=j.as1[0][5]; a32=j.aa3[26]*j.cep6; a12=j.aa1[26]*j.cmp6; s12=j.sa1[26]*j.csp6; break;
    	case "D33(1700)": am=j.am3[1];    ae=j.ae3[1];    as=j.as3[1];    a32=j.aa3[ 2]*j.cep7; a12=j.aa1[ 2]*j.cmp7; s12=j.sa1[ 2]*j.csp7; break;
    	case "P13(1720)": am=j.am1[0][6]; ae=j.ae1[0][6]; as=j.as1[0][6]; a32=j.aa3[27]*j.cep8; a12=j.aa1[27]*j.cmp8; s12=j.sa1[27]*j.csp8; break;
    	case "D15(1675)": am=j.am1[0][8]; ae=j.ae1[0][8]; as=0;           a32=j.aa3[29]*j.cep9; a12=j.aa1[29]*j.cmp9; s12=0;                break;
    	case "F35(1905)": am=j.am3[5];    ae=j.ae3[5];    as=0;           a32=j.aa3[ 6]*j.cep10;a12=j.aa1[ 6]*j.cmp10;s12=0;                break;
    	case "F37(1950)": am=j.am3[2];    ae=j.ae3[2];    as=0;           a32=j.aa3[ 3]*j.cep11;a12=j.aa1[ 3]*j.cmp11;s12=0; 
    	}
    	double[] dcof = {am,ae,as,a32,a12,s12}; String[] scof = new String[6];
    	int n=0; for (double d : dcof) scof[n++] = String.valueOf(""+String.format("%.2f", d));
    	setCofParameters(scof);
    }
    
    public void setCofParameters(String[] val) {
        cof_parameters.setText("<html>"+
                "<b>MULTIPOLE:</b>&nbsp;AM="+val[0]+"&nbsp;&nbsp;&nbsp;&nbsp;AE="+val[1]+"&nbsp;&nbsp;AS="+val[2]+"<br>" +
                "<b>ISOSPIN:</b>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;A32="+val[3]+"&nbsp;A12="+val[4]+"&nbsp;S12="+val[5]+"<br>" +
                "</html>"
            );
    }
    
    public JButton getResetButton() {
        JButton button = new JButton("Reset");
        button.addActionListener(this); 
        return button;
    }

    public void setJanrTabNames(int i, String... names) {
        for(String name : names) detectorTabNames[i].add(name); 
        if(dgmActive) {dgm[i].setDetectorTabNames(names); setDetectorCanvas(i,dgm[i].getDetectorCanvas()); return;}
        setDetectorCanvas(i,new EmbeddedCanvasTabbed(names));      
    }
    
    public void setDetectorCanvas(int i, EmbeddedCanvasTabbed canvas) {
        detectorCanvas[i] = canvas;
    }
    
    public void setDetectorCanvas(EmbeddedCanvasTabbed canvas) {
        detectorCanvas[0] = canvas;
    }
    
    public EmbeddedCanvasTabbed getDetectorCanvas(int i) {
        return detectorCanvas[i];
    }    
    
    public JPanel getDetectorPanel() {
        return detectorPanel;
    }
    
    public String getDetectorName() {
        return detectorName;
    }
    
    public void createHistos(int run) {

    }
    
    public void plotHistos(int run, int nt, int bin) {
  	
    } 
    
    public void hjanr_loadpar(int n) {
    	
    }
    
    public int getRunNumber() {
        return runNumber;
    }
    
    public void setRunNumber(int val) {
    	runNumber = val;
    }
    
    public void timerUpdate() { 
    	
    }
    
    public void ResetAction() {
    	cofList.clear(); configCMB.setSelectedIndex(0);
    	cofslider.setValue(((int)j.start_value[0]+5)*10);  
    	hjanr_loadpar(3); 
    	plotHistos(getRunNumber(),nt,bin); 
    	Res2Amp(j.res_name[cof]);    	
    }
    
    public void resetEventListener() {
        System.out.println(root+"resetEventListener():" +  getDetectorName() + " histogram for run "+ getRunNumber());
        createHistos(100);
        plotHistos(100,-1,-1);
    }

    public void setCanvasUpdate(int time) {    	
    	for(int i=0; i<2; i++) {
    		for(int tab=0; tab<detectorTabNames[i].size(); tab++) {
    			System.out.println(detectorTabNames[i].get(tab));
    			getDetectorCanvas(i).getCanvas(detectorTabNames[i].get(tab)).initTimer(time);
    		}
    	}
    }

    public void actionPerformed(ActionEvent e) { 
    	if(e.getActionCommand().compareTo("Reset")==0) ResetAction();
    	if(bG0!=null)  binWCF = Integer.parseInt(bG0.getSelection().getActionCommand()); 
        if(bG1!=null) bin1234 = Integer.parseInt(bG1.getSelection().getActionCommand());
        nt = bin1234 + 4*binWCF;
        plotHistos(getRunNumber(),nt,bin);
    }
}

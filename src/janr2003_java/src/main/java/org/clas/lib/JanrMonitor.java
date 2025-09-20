package org.clas.lib;

import java.awt.BorderLayout;
import java.awt.Dimension;
import java.awt.FlowLayout;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.ItemEvent;
import java.awt.event.ItemListener;
import java.util.ArrayList;

import javax.swing.BorderFactory;
import javax.swing.ButtonGroup;
import javax.swing.JButton;
import javax.swing.JCheckBox;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JRadioButton;
import javax.swing.JSlider;
import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;

import org.jlab.groot.base.GStyle;

public class JanrMonitor implements ActionListener {
	
    private String                 detectorName = null;
    private JPanel                detectorPanel = null;
    private EmbeddedCanvasTabbed detectorCanvas = null;
    private ArrayList<String>  detectorTabNames = new ArrayList();
    public  DataGroupManager                dgm = new DataGroupManager();
    public Boolean                    dgmActive = false;
    public int                        runNumber = 0;

    
    private int bin1234                         = 0;
    private int binWCF                          = 0;
    public  int nt                              = 0;
    private int bin                             = 0;
    		
    private JPanel                  actionPanel = null;
    private JPanel               controlsPanel0 = null;
    private JPanel               controlsPanel1 = null;
 
    public Boolean                  histosExist = false;   
    public String                          root = " ";
    
    public  JSlider                   binslider = null;
    
    private ButtonGroup bG0, bG1;
    private JRadioButton bW,bC,bF,b1,b2,b3,b4;
    private Boolean useWCF=false, use1234=false, useBSP=false;
    
    public JanrMonitor(String name) {
		initGStyle(14);
	    detectorName   = name;
	    detectorPanel  = new JPanel();
	    detectorCanvas = new EmbeddedCanvasTabbed();
        root           = detectorName+".JanrMonitor.";		
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
    
    public void initPanel() {
        getDetectorPanel().setLayout(new BorderLayout());
        actionPanel    = new JPanel(new FlowLayout());
        controlsPanel0 = new JPanel(new GridBagLayout());
        
        controlsPanel1 = new JPanel();
        controlsPanel1.setBorder(BorderFactory.createTitledBorder("Display"));		       
        controlsPanel1.add(packActionPanel());
        controlsPanel1.add(packBinSliderPanel());

        GridBagConstraints c = new GridBagConstraints();
        c.fill = GridBagConstraints.HORIZONTAL; c.weightx = 0.5;		
        c.gridx=0 ; c.gridy=0 ; controlsPanel0.add(controlsPanel1,c);

        getDetectorPanel().add(getDetectorCanvas(),BorderLayout.CENTER);           
        getDetectorPanel().add(controlsPanel0,     BorderLayout.SOUTH); 
    }
    
    public JPanel packActionPanel() {
    	if (useWCF || use1234) actionPanel.add(getButtonPanel());
    	return actionPanel;
    }
    
    public JPanel packBinSliderPanel() {
        actionPanel.add(getBinSliderPane());	
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
    
    public JPanel getButtonPanel() {
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
    
    public JPanel getBinSliderPane() {
        JPanel sliderPane = new JPanel();
        JLabel      label = new JLabel("" + String.format("%d", 0));       
        binslider         = new JSlider(JSlider.HORIZONTAL, 0, 30, 10); 
        binslider.setPreferredSize(new Dimension(300,10));
        sliderPane.add(new JLabel("WCF bin",JLabel.CENTER));
        sliderPane.add(binslider);
        sliderPane.add(label);  
        binslider.addChangeListener(new ChangeListener() {
            public void stateChanged(ChangeEvent e) {
                JSlider slider = (JSlider) e.getSource(); 
                bin = slider.getValue(); 
                label.setText(String.valueOf(""+String.format("%d", bin)));
//                System.out.println(detectorCanvas.getCanvas().getClickedPad()+" "+binWCF+" "+bin1234+" "+bin);
//                int nt = detectorCanvas.getCanvas().getClickedPad();
                plotHistos(getRunNumber(), nt, bin);
            }
        });        
        return sliderPane;
    }    
    
    public void setJanrTabNames(String... names) {
        for(String name : names) detectorTabNames.add(name); 
        if(dgmActive) {dgm.setDetectorTabNames(names); setDetectorCanvas(dgm.getDetectorCanvas()); return;}
        setDetectorCanvas(new EmbeddedCanvasTabbed(names));      
    }
    
    public void setDetectorCanvas(EmbeddedCanvasTabbed canvas) {
        detectorCanvas = canvas;
    }
    
    public EmbeddedCanvasTabbed getDetectorCanvas() {
        return detectorCanvas;
   }
    
    public JPanel getDetectorPanel() {
        return detectorPanel;
    }
    
    public String getDetectorName() {
        return detectorName;
    }
    
    public void createHistos(int run) {
        // initialize canvas and create histograms
    }
    
    public void plotHistos(int run, int nt, int bin) {
  	
    } 
    
    public int getRunNumber() {
        return runNumber;
    }
    
    public void setRunNumber(int val) {
    	runNumber = val;
    }
    
    public void timerUpdate() { 
    	
    }
    
    public void resetEventListener() {
        System.out.println(root+"resetEventListener():" +  getDetectorName() + " histogram for run "+ getRunNumber());
        createHistos(100);
        plotHistos(100,-1,-1);
    }
    
    public void setCanvasUpdate(int time) {
        for(int tab=0; tab<detectorTabNames.size(); tab++) {
           getDetectorCanvas().getCanvas(detectorTabNames.get(tab)).initTimer(time);
        }
    }
    
    public void actionPerformed(ActionEvent e) {
        
    	if(bG0!=null)  binWCF = Integer.parseInt(bG0.getSelection().getActionCommand()); 
        if(bG1!=null) bin1234 = Integer.parseInt(bG1.getSelection().getActionCommand());
        nt = bin1234 + 4*binWCF;
//        System.out.println(nt+" "+bin);
        plotHistos(getRunNumber(),nt,bin);
    }
}

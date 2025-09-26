package org.clas.lib;

import java.awt.BorderLayout;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.ItemEvent;
import java.awt.event.ItemListener;

import javax.swing.JFileChooser;
import javax.swing.JFrame;
import javax.swing.JMenu;
import javax.swing.JMenuBar;
import javax.swing.JMenuItem;
import javax.swing.JPanel;
import javax.swing.JTabbedPane;
import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;
import javax.swing.filechooser.FileSystemView;

import org.jlab.detector.view.DetectorListener;
import org.jlab.detector.view.DetectorShape2D;
import org.jlab.groot.base.GStyle;
import org.jlab.io.base.DataEvent;
import org.jlab.io.task.IDataEventListener;

public class JanrViewer implements IDataEventListener, DetectorListener, ActionListener, ItemListener, ChangeListener {

    JFrame frame = new JFrame("JANR2003");

	JPanel                   mainPanel = new JPanel();	
	JMenuBar                   menuBar = new JMenuBar();     
    JFileChooser                    fc = null; 	
	JanrMonitor[]             monitors = null;

	java.util.Timer       processTimer = null;
    int                     eventDelay = 0;
    int               canvasUpdateTime = 2000;
    int             analysisUpdateEvnt = 100;
    int                      runNumber = 1;;

    String outPath = "/Users/cole/JANR2025/";
    String workDir = outPath;
    	
    int    selectedTabIndex = 0;
    String selectedTabName  = " ";
    	    
	public static void main(String[] args){
		System.out.println("*** WELCOME TO JANR2003 ***\n");
	    JanrViewer viewer = new JanrViewer(args);
	}
	
    public JanrViewer(String[] args) {   
    	createMonitors(args);
    	createMenuBar();
    	createPanels();  
	    frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
	    frame.setJMenuBar(getMenuBar());
	    frame.add(getMainPanel());
	    frame.pack();
	    frame.setSize(2300, 1100);
	    frame.setVisible(true);
        System.out.println("JanrViewer init complete \n");
    }

    public void createMonitors(String[] args) {
        workDir = FileSystemView.getFileSystemView().getHomeDirectory().toString()+"/ECAL/";
    	System.out.println("JanrViewer.createMonitors: workDir="+workDir);
     	monitors = new JanrMonitor[args.length==0 ? 1:args.length];
        monitors[0] = new HJanr("JANR");
        initMonitors();
    }
    
    public void initMonitors() {
    	System.out.println("JanrViewer.initMonitors");
//		for(int k=0; k<this.monitors.length; k++) monitors[k].init();
    } 
    
    public void createMenuBar() {
    	System.out.println("JanrViewer.createMenuBar");        		         
        JMenu menu = new JMenu("File");
        JMenuItem menuItem = new JMenuItem("Fit Data"); menuItem.addActionListener(this); menu.add(menuItem);
        getMenuBar().add(menu);
    }
    
    public JMenuBar getMenuBar() {
    	return menuBar;
    }
    
    public JPanel getMainPanel() {
    	return mainPanel;
    }
    
    public void createPanels() {
    	System.out.println("JanrViewer.createPanels()");
    	
      	JTabbedPane tabbedpane = new JTabbedPane(); tabbedpane.addChangeListener(this);	               
        getMainPanel().setLayout(new BorderLayout()); getMainPanel().add(tabbedpane);
       
        GStyle.getAxisAttributesX().setTitleFontSize(18);
        GStyle.getAxisAttributesX().setLabelFontSize(14);
        GStyle.getAxisAttributesY().setTitleFontSize(18);
        GStyle.getAxisAttributesY().setLabelFontSize(14);

        for(int k=0; k<this.monitors.length; k++) tabbedpane.add(this.monitors[k].getDetectorPanel(), this.monitors[k].getDetectorName());                     

        setCanvasUpdate(canvasUpdateTime);       
    }
        
    public void setCanvasUpdate(int time) {
        System.out.println("JanrViewer.setCanvasUpdate("+time+")");
        this.canvasUpdateTime = time;
        for(int k=0; k<this.monitors.length; k++) {
            this.monitors[k].setCanvasUpdate(time);
        }
    }
    
	@Override
	public void itemStateChanged(ItemEvent e) {
		
	}
    
    @Override
    public void actionPerformed(ActionEvent e) {

    }
    
    @Override
    public void timerUpdate() {
        if(this.runNumber==0) return;
        for(int k=0; k<this.monitors.length; k++) {
            this.monitors[k].timerUpdate();
        }
   }
   
    @Override
    public void dataEventAction(DataEvent event) {  
    	
    }
    
    @Override
    public void resetEventListener() {
        for(int k=0; k<this.monitors.length; k++) {
            this.monitors[k].resetEventListener();
            this.monitors[k].timerUpdate();
        }      
    } 
   
	@Override
	public void processShape(DetectorShape2D arg0) {
		// TODO Auto-generated method stub		
	}
    
    @Override
    public void stateChanged(ChangeEvent e) {
        this.timerUpdate();
    }    
    
}

package org.vardb.util;

import java.awt.BorderLayout;
import java.util.Random;

import javax.swing.JFrame;

import org.vardb.util.CException;
import org.vardb.util.DataFrame;

import weka.classifiers.Classifier;
import weka.classifiers.Evaluation;
import weka.classifiers.evaluation.ThresholdCurve;
import weka.classifiers.functions.Logistic;
import weka.classifiers.trees.J48;
import weka.core.Attribute;
import weka.core.FastVector;
import weka.core.Instance;
import weka.core.Instances;
import weka.core.Utils;
import weka.core.converters.ConverterUtils.DataSink;
import weka.core.converters.ConverterUtils.DataSource;
import weka.gui.visualize.PlotData2D;
import weka.gui.visualize.ThresholdVisualizePanel;

public class WekaHelper
{
	 public static void test(String filename)
	 {
		 try
		 {
			 Instances data = loadData(filename);
			 System.out.println(data.toSummaryString());
			 //J48 classifier=buildJ48Classifier(data);
			 Classifier classifier=buildLogisticClassifier(data);
			 Evaluation eval=crossValidate(data,classifier,10);
			 getReport(eval);
			 plotRocCurve(eval);
		 }
		 catch(Exception e)
		 {
			 throw new CException(e);
		 }
	 }
	 
	 public static void getReport(Evaluation eval)
	 {
		 try
		 {
			System.out.println(eval.toSummaryString() + "\n");
			System.out.println(eval.toMatrixString() + "\n");
			System.out.println(eval.toClassDetailsString() + "\n");
		 }
		 catch(Exception e)
		 {
			 throw new CException(e);
		 }
	 }
	 
	 public static Evaluation crossValidate(Instances data, Classifier classifier, int folds)
	 {
		 try
		 {
			 Evaluation eval = new Evaluation(data);
			 eval.crossValidateModel(classifier, data, folds, new Random(1));
			 return eval;
		 }
		 catch(Exception e)
		 {
			 throw new CException(e);
		 }
	 }
	 
	 public static Logistic buildLogisticClassifier(Instances data)
	 {
		 try
		 {
			 Logistic classifier=new Logistic();
			 classifier.buildClassifier(data);
			 System.out.println(classifier.toString());
			 System.out.println("coefficients: "+classifier.coefficients());
			 return classifier;
		 }
		 catch(Exception e)
		 {
			 throw new CException(e);
		 }
	 } 
	 
	 public static J48 buildJ48Classifier(Instances data)
	 {
		 try
		 {
			 J48 tree = new J48();
			 tree.setUnpruned(false);
			 tree.buildClassifier(data);
			 System.out.println(tree.toString());
			 //System.out.println(tree.toSummaryString() + "\n");
			 return tree;
		 }
		 catch(Exception e)
		 {
			 throw new CException(e);
		 }
	 }
	 
	 public static Instances loadData(String filename)
	 {
		 try
		 {
			 DataSource source = new DataSource(filename);
			 Instances data = source.getDataSet();
			 // setting class attribute if the data format does not provide this information
			 // E.g., the XRFF format saves the class attribute information as well
			 if (data.classIndex() == -1)
			   data.setClassIndex(data.numAttributes() - 1);
			 return data;
			 
		 }
		 catch(Exception e)
		 {
			 throw new CException(e);
		 }
	 }
	 
	 public static Attribute createIntAttribute(String name)
	 {
		 return new Attribute(name);
	 }
	 
	 public static Attribute createNominalAttribute(String name, String...args)
	 {
		 FastVector labels = new FastVector();
		 for (String arg : args)
		 {
			 labels.addElement(arg);
		 }
		 return new Attribute(name,labels);
	 }
	 
	 public static Instances convertInstances(DataFrame<?,?,?> dataframe)
	 {
		 FastVector attributes = new FastVector();
		 for (Object columnKey : dataframe.columnKeySet())
		 {
			 Attribute attribute=createIntAttribute(columnKey.toString());
			 attributes.addElement(attribute);
		 }
		 Instances data = new Instances("dataset", attributes, dataframe.getNumRows());
		 for (Object rowKey : dataframe.rowKeySet())
		 {
			 double[] values = new double[data.numAttributes()];
			 int index=0;
			 for (Object columnKey : dataframe.columnKeySet())
			 {
				 Object cell = dataframe.get(rowKey, columnKey);
				 values[index]=Double.parseDouble(cell.toString());
				 Instance instance = new Instance(1.0, values);
				 data.add(instance);
				 index++;
			 }
		 }
		 return data;
	 }
	 
	 public static void writeData(Instances data, String filename)
	 {
		try
		{
			 DataSink.write(filename, data);
		}
		catch(Exception e)
		{
			throw new CException(e);
		}
	 }
	 
	 public static void plotRocCurve(Evaluation eval)
	 {
		try
		{
			ThresholdCurve tc = new ThresholdCurve();
			int classIndex = 0; // ROC for the 1st class label
			Instances curve = tc.getCurve(eval.predictions(), classIndex);
			
			//2. Put the plotable into a plot container
			PlotData2D plotdata = new PlotData2D(curve);
			plotdata.setPlotName(curve.relationName());
			plotdata.addInstanceNumberAttribute();
			
			//3. Add the plot container to a visualization panel
			ThresholdVisualizePanel tvp = new ThresholdVisualizePanel();
			tvp.setROCString("(Area under ROC = " +
					Utils.doubleToString(ThresholdCurve.getROCArea(curve),4)+")");
			tvp.setName(curve.relationName());
			tvp.addPlot(plotdata);
			
			//4. Add the visualization panel to a JFrame
			final JFrame jf = new JFrame("WEKA ROC: " + tvp.getName());
			jf.setSize(500,400);
			jf.getContentPane().setLayout(new BorderLayout());
			jf.getContentPane().add(tvp, BorderLayout.CENTER);
			jf.setDefaultCloseOperation(JFrame.DISPOSE_ON_CLOSE);
			jf.setVisible(true);
		}
		catch(Exception e)
		{
			throw new CException(e);
		}
	 }
}
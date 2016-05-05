package com.company;
import com.company.Point;
import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartPanel;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.data.category.DefaultCategoryDataset;
import org.jfree.ui.ApplicationFrame;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;
import org.jfree.ui.RefineryUtilities;

import java.util.ArrayList;
import java.util.Vector;

/**
 * Created by admin on 01/10/15.
 */
public class Plotter extends ApplicationFrame {
    public Plotter(final String title, ArrayList<Vector<Point>> data, String [] names) {

        super(title);

        final XYSeriesCollection data_plot = new XYSeriesCollection();
        for (int j=0; j<data.size(); j++) {
            final XYSeries series = new XYSeries(names[j]);
            for (int i = 0; i < data.get(j).size(); i++) {
                series.add(data.get(j).get(i).x, data.get(j).get(i).y);
            }

            data_plot.addSeries(series);
        }
        final JFreeChart chart = ChartFactory.createXYLineChart(
                "",
                "X",
                "Y",
                data_plot,
                PlotOrientation.VERTICAL,
                true,
                true,
                false
        );

        final ChartPanel chartPanel = new ChartPanel(chart);
        chartPanel.setPreferredSize(new java.awt.Dimension(500, 270));
        setContentPane(chartPanel);

    }

    void plot(){
        pack();
        RefineryUtilities.centerFrameOnScreen(this);
        setVisible(true);
    }
}

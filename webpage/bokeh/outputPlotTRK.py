#----- Author: Nick Konz -----#

from bokeh.layouts import column, row, widgetbox, gridplot
from bokeh.models import Band, CustomJS, ColumnDataSource, Dropdown, Slider, Button, BoxZoomTool, Range1d, BoxSelectTool, LassoSelectTool, Select, CheckboxGroup, Whisker
from bokeh.plotting import Figure, output_file, show, save, curdoc
from bokeh.transform import linear_cmap
from bokeh.models.tools import HoverTool
from bokeh.palettes import Viridis256
from bokeh.events import MouseEnter

import numpy as np

output_file("outputPlotTRK.html")

x = []
y = []
err_x_up = []
err_y_up = []
err_x_down = []
err_y_down = []

x_m = []
y_m = []
lower1 = []
lower2 = []
lower3 = []
upper1 = []
upper2 = []
upper3 = []


# TESTING ################################################

#x = [1,4]
#y = [-1,4]

#err_x = [2.0, 0.3] #first is right x EB on left point; second is same but on right point.
#err_y = [3.4, 1.6] #first is BOTH y EB on left point; second is same but on right point.

#err_x_up = [x[0] + err_x[0], x[1] + err_x[1]]
#err_y_up = [y[0] + err_y[0], y[1] + err_y[1]]
#err_x_down = [x[0] - err_x[0], x[1] - err_x[1]]
#err_y_down = [y[0] - err_y[0], y[1] - err_y[1]]

#m = 5/3
#b = -8/3
#def testlin(x):
#    return m*x + b


#N = 1000

#x_m = np.linspace(min(np.array(err_x_down)) - 1, max(np.array(err_x_up )) + 1, N)
#y_m = testlin(x_m)

#x_m = x_m.tolist()
#y_m = y_m.tolist()

#testSlopy = 1.0

#x_m = np.array(x_m)
#y_m = np.array(y_m)
#lower1 = y_m - testSlopy
#upper1 = y_m + testSlopy
#lower2 = y_m - 2*testSlopy
#upper2 = y_m + 2*testSlopy
#lower3 = y_m - 3*testSlopy
#upper3 = y_m + 3*testSlopy

#x_m.tolist()
#y_m.tolist()
#lower1.tolist()
#upper1.tolist()
#lower2.tolist()
#upper2.tolist()
#lower3.tolist()
#upper3.tolist()

# ################################################

srcData = ColumnDataSource(data=dict(x = x, y = y, 
                                     err_x_up = err_x_up, 
                                     err_y_up = err_y_up, 
                                     err_x_down = err_x_down, 
                                     err_y_down = err_y_down
                                     ))
srcModel = ColumnDataSource(data=dict(x_m = x_m, y_m = y_m, 
                                      lower1 = lower1, 
                                      lower2 = lower2,
                                      lower3 = lower3,
                                      upper1 = upper1,
                                      upper2 = upper2,
                                      upper3 = upper3
                                      ))

# create the scatter plot
TOOLS="pan,wheel_zoom,reset"

p = Figure(tools=TOOLS, plot_width=600, plot_height=600, min_border=10, min_border_left=50,
           title="Fitted Model",
           toolbar_location="above",
           active_scroll='wheel_zoom', 
           active_drag = "pan",
           x_axis_label = "x",
           y_axis_label = "y"
           )
p.background_fill_color = "#fafafa"

r = p.scatter(source = srcData, x='x', y='y', size=10, color="red", alpha=1.0)

#y error bars:
p.add_layout(
        Whisker(source = srcData, base="x", upper="err_y_up", lower="err_y_down")
    )

#x error bars:
p.add_layout(
        Whisker(source = srcData, base="y", upper="err_x_up", lower="err_x_down", dimension="width")
    )

#model plot

p.line(source=srcModel, x = "x_m", y = "y_m", line_width=3, line_alpha=0.6, color="black")

#slop envelopes
p.add_layout(
        Band(source = srcModel, base = "x_m", lower="lower1", upper="upper1", level="underlay",
             fill_color = "blue", fill_alpha=0.3, line_width=2, line_color="black")
    )
p.add_layout(
        Band(source = srcModel, base = "x_m", lower="lower2", upper="upper2", level="underlay",
             fill_color = "blue", fill_alpha=0.2, line_width=2, line_color="black")
    )
p.add_layout(
        Band(source = srcModel, base = "x_m", lower="lower3", upper="upper3", level="underlay",
             fill_color = "blue", fill_alpha=0.1, line_width=2, line_color="black")
    )

#JS callbacks

callbackPlot = CustomJS(args=dict(srcData=srcData, srcModel=srcModel, p=p, xaxis=p.xaxis[0], yaxis=p.yaxis[0]), code="""
    p.reset.emit();

    let input_data = getData();
    let model_plot = getModelData();

    let x_data = input_data[0];
    let err_x_data = input_data[1];
    let y_data = input_data[2];
    let err_y_data = input_data[3];

    let x_m_data = model_plot[0];
    let y_m_data = model_plot[1];
    let prams = model_plot[3];

    let x = [];
    let err_x_up = [];
    let err_x_down = [];
    let y = []; 
    let err_y_up = [];
    let err_y_down = [];

    let x_m = [];
    let y_m = [];
    let lower1 = [];
    let lower2 = [];
    let lower3 = [];
    let upper1 = [];
    let upper2 = [];
    let upper3 = [];
    
    var dataData = srcData.data;
    var dataModel = srcModel.data;

    dataData['x'] = [];
    dataData['y'] = [];
    dataData['err_x_up'] = [];
    dataData['err_y_up'] = [];
    dataData['err_x_down'] = [];
    dataData['err_y_down'] = [];

    dataModel['x_m'] = [];
    dataModel['y_m'] = [];
    dataModel['lower1'] = [];
    dataModel['lower2'] = [];
    dataModel['lower3'] = [];
    dataModel['upper1'] = [];
    dataModel['upper2'] = [];
    dataModel['upper3'] = [];

    x = dataData['x'];
    y = dataData['y'];
    err_x_up = dataData['err_x_up']
    err_y_up = dataData['err_y_up']
    err_x_down = dataData['err_x_down']
    err_y_down = dataData['err_y_down']
    
    x_m = dataModel['x_m'];
    y_m = dataModel['y_m'];
    lower1 = dataModel['lower1'];
    lower2 = dataModel['lower2'];
    lower3 = dataModel['lower3'];
    upper1 = dataModel['upper1'];
    upper2 = dataModel['upper2'];
    upper3 = dataModel['upper3'];

    let bands = getBands(y_m_data, prams);

    for (let i = 0; i < x_data.length; i++){
        x[i] = x_data[i];
        y[i] = y_data[i];
        err_x_up[i] = x_data[i] + err_x_data[i];
        err_y_up[i] = y_data[i] + err_y_data[i];
        err_x_down[i] = x_data[i] - err_x_data[i];
        err_y_down[i] = y_data[i] - err_y_data[i];
    }

    for (let i = 0; i < x_m_data.length; i++){
        x_m[i] = x_m_data[i];
        y_m[i] = y_m_data[i];

        lower1[i] = bands[0][i];
        lower2[i] = bands[1][i];
        lower3[i] = bands[2][i];
        upper1[i] = bands[3][i];
        upper2[i] = bands[4][i];
        upper3[i] = bands[5][i];
    }

    dataData['x'] = x;
    dataData['y'] = y;
    dataData['err_x_up'] = err_x_up;
    dataData['err_y_up'] = err_y_up;
    dataData['err_x_down'] = err_x_down;
    dataData['err_y_down'] = err_y_down;

    dataModel['x_m'] = x_m;
    dataModel['y_m'] = y_m;
    dataModel['lower1'] = lower1;
    dataModel['lower2'] = lower2;
    dataModel['lower3'] = lower3;
    dataModel['upper1'] = upper1;
    dataModel['upper2'] = upper2;
    dataModel['upper3'] = upper3;

    console.log(lower1)
    console.log(lower2)
    console.log(lower3)
    console.log(upper1)
    console.log(upper2)
    console.log(upper3)

    srcData.change.emit();  
    srcModel.change.emit(); 
    p.change.emit();
""")


layout = p

#interactivity

layout.js_on_event(MouseEnter, callbackPlot) #plot whichever axes are currently selected

curdoc().add_root(layout) #any changes to layout will trigger on_change callbacks
curdoc().title = "Fitted Model"

save(layout)
#show(layout)
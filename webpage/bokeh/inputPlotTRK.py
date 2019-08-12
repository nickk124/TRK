#----- Author: Nick Konz -----#

from bokeh.layouts import column, row, widgetbox, gridplot
from bokeh.models import CustomJS, ColumnDataSource, Dropdown, Slider, Button, BoxZoomTool, Range1d, BoxSelectTool, LassoSelectTool, Select, CheckboxGroup, Whisker
from bokeh.plotting import Figure, output_file, show, save, curdoc
from bokeh.transform import linear_cmap
from bokeh.models.tools import HoverTool
from bokeh.palettes import Viridis256
from bokeh.events import MouseEnter

import numpy as np

output_file("inputPlotTRK.html")

x = []
y = []
err_x_up = []
err_y_up = []
err_x_down = []
err_y_down = []

#x = [1,4]
#y = [-1,4]

#err_x = [2.0, 0.3] #first is right x EB on left point; second is same but on right point.
#err_y = [3.4, 1.6] #first is BOTH y EB on left point; second is same but on right point.

#err_x_up = [x[0] + err_x[0], x[1] + err_x[1]]
#err_y_up = [y[0] + err_y[0], y[1] + err_y[1]]
#err_x_down = [x[0] - err_x[0], x[1] - err_x[1]]
#err_y_down = [y[0] - err_y[0], y[1] - err_y[1]]

srcData = ColumnDataSource(data=dict(x = x, y = y, err_x_up = err_x_up, err_y_up = err_y_up, err_x_down = err_x_down, err_y_down = err_y_down))


# create the scatter plot
TOOLS="pan,wheel_zoom,reset"

p = Figure(tools=TOOLS, plot_width=500, plot_height=500, min_border=10, min_border_left=50,
           toolbar_location="above",
           active_scroll='wheel_zoom', 
           active_drag = "pan"
           )
p.background_fill_color = "#fafafa"

r = p.scatter(source = srcData, x='x', y='y', size=10, color="blue", alpha=0.6)

#y error bars:
p.add_layout(
        Whisker(source = srcData, base="x", upper="err_y_up", lower="err_y_down")
    )

#x error bars:
p.add_layout(
        Whisker(source = srcData, base="y", upper="err_x_up", lower="err_x_down", dimension="width")
    )

#p.add_tools(HoverTool(tooltips=[("x", "@x"), ("err_x", "@err_x"), ("y", "@y"), ("err_y", "@err_y")]))

#JS callbacks

callbackPlot = CustomJS(args=dict(srcData=srcData, p=p, xaxis=p.xaxis[0], yaxis=p.yaxis[0]), code="""
    p.reset.emit();

    let input_data = getData();

    let x_data = input_data[0];
    let err_x_data = input_data[1];
    let y_data = input_data[2];
    let err_y_data = input_data[3];

    let x = [];
    let err_x_up = [];
    let err_x_down = [];
    let y = []; 
    let err_y_up = [];
    let err_y_down = [];
    
    var dataData = srcData.data;

    dataData['x'] = [];
    dataData['y'] = [];
    dataData['err_x_up'] = [];
    dataData['err_y_up'] = [];
    dataData['err_x_down'] = [];
    dataData['err_y_down'] = [];

    x = dataData['x'];
    y = dataData['y'];
    err_x_up = dataData['err_x_up']
    err_y_up = dataData['err_y_up']
    err_x_down = dataData['err_x_down']
    err_y_down = dataData['err_y_down']

    for (let i = 0; i < x_data.length; i++){
        x[i] = x_data[i];
        y[i] = y_data[i];
        err_x_up[i] = x_data[i] + err_x_data[i];
        err_y_up[i] = y_data[i] + err_y_data[i];
        err_x_down[i] = x_data[i] - err_x_data[i];
        err_y_down[i] = y_data[i] - err_y_data[i];
    }

    dataData['x'] = x;
    dataData['y'] = y;
    dataData['err_x_up'] = err_x_up;
    dataData['err_y_up'] = err_y_up;
    dataData['err_x_down'] = err_x_down;
    dataData['err_y_down'] = err_y_down;

    //xaxis.axis_label = playerPlotAxes[0];
    //yaxis.axis_label = playerPlotAxes[1];

    srcData.change.emit();   
    p.change.emit();
""")

#callbackUpdateHAxis = CustomJS(code="""
#    playerPlotAxes[0] = this.value;
#""")

#callbackUpdateVAxis = CustomJS(code="""
#    playerPlotAxes[1] = this.value;
#""")

#widgets


#hOptions = ["KO Counts per Player", "Win Counts per Player", "Damage Done per Player"]
#vOptions = ["KO Counts per Player", "Win Counts per Player", "Damage Done per Player"]

#hSelect = Select(title =    "Horizontal Data", 
#                 options =  hOptions,   #list of options
#                 value =   "KO Counts per Player" #default value  
#                 )

#vSelect = Select(title =    "Vertical Data", 
#                 options =  vOptions,   #list of options
#                 value =   "Win Counts per Player" #default value  
#                 )

#layout

#widgets = row(hSelect, vSelect)
#layout = column(widgets, p)
layout = p

#interactivity

layout.js_on_event(MouseEnter, callbackPlot) #plot whichever axes are currently selected

#hSelect.js_on_change('value', callbackUpdateHAxis)
#vSelect.js_on_change('value', callbackUpdateVAxis)

#hSelect.js_on_change('value', callbackPlot)
#vSelect.js_on_change('value', callbackPlot)

curdoc().add_root(layout) #any changes to layout will trigger on_change callbacks
curdoc().title = "Input Data"

save(layout)
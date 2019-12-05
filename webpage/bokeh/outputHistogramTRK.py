from bokeh.layouts import column, row, widgetbox, gridplot
from bokeh.models import Band, CustomJS, ColumnDataSource, Dropdown, Slider, Button, BoxZoomTool, Range1d, BoxSelectTool, LassoSelectTool, Select, CheckboxGroup, Whisker
from bokeh.plotting import Figure, output_file, show, save, curdoc
from bokeh.transform import linear_cmap
from bokeh.models.tools import HoverTool
from bokeh.palettes import Viridis256
from bokeh.events import MouseEnter

import numpy as np

output_file("outputHistogramTRK.html")

"""
y = [1,1,1,2,3,4,9,10,11,15,16,17,18,55]

N = len(y)

bincount = np.int(np.sqrt(float(N)))
range = [np.min(y), np.max(y)]

arr_hist, edges = np.histogram(y, bins = bincount, range = range ) #bins is the number of bins
left = edges[:-1]
right = edges[1:]
print(edges)
print(bincount)
print(range[1] - range[0])
"""

arr_hist = []
left = []
right = []
xMin = 0
xMax = 1

src = ColumnDataSource(data=dict(arr_hist = arr_hist, left = left, right = right))

p = Figure(plot_height = 400, plot_width = 500, x_range = (xMin, xMax),
                    x_axis_label = 'parameter',
                    y_axis_label = 'probability', tools = "xpan, xwheel_zoom", active_scroll='xwheel_zoom', active_drag = "xpan")

p.quad(source = src, bottom = 0, top = 'arr_hist', left = 'left', right = 'right', fill_color = 'cornflowerblue', line_color = 'black')

#defines the callback to be used:
callback_plot = CustomJS(args=dict(src=src, p=p, axis=p.xaxis[0], x_range=p.x_range), code="""

    p.reset.emit();


    var hist_result = paramHistograms[whichParam];

    var arr_hist_result = hist_result[0];
    var left_result = hist_result[1].slice(0,-1);
    var right_result = hist_result[1].slice(1);

    let bincount = arr_hist_result.length;

    var data = src.data;

    var arr_hist = [];
    var left = [];
    var right = [];

    data['arr_hist'] = [];
    data['left'] = [];
    data['right'] = [];

    arr_hist = data['arr_hist'];
    left = data['left'];
    right = data['right'];

    for (var i = 0; i < bincount; i++){
        arr_hist[i] = arr_hist_result[i];
        left[i] = left_result[i];
        right[i] = right_result[i];
    }

    data['arr_hist'] = arr_hist;
    data['left'] = left;
    data['right'] = right;

    x_range.start = left[0];
    x_range.end = right[bincount-1];

    p.change.emit();
    src.change.emit();

""")

#callback_addbins = CustomJS(args=dict(src=src, p=p, x_range=p.x_range), code="""

#    var oldxMin = xMin;
#    var oldxMax = xMax;

#    p.reset.emit();

#    hasWeights = $('#Weighted').data('clicked');

#    var data_res = getData();
#    var y_data = data_res[0];
#    var w_data = [];

#    if (hasWeights){
#        w_data = data_res[1];
#    }


#    bincount += 1;

#    var hist_result = getHistBins(y_data, w_data);
#    var arr_hist_result = hist_result[0];
#    var left_result = hist_result[1];
#    var right_result = hist_result[2];

#    var data = src.data;

#    var arr_hist = [];
#    var left = [];
#    var right = [];

#    arr_hist = data['arr_hist'];
#    left = data['left'];
#    right = data['right'];

#    for (var i = 0; i < bincount; i++){
#        arr_hist[i] = arr_hist_result[i];
#        left[i] = left_result[i];
#        right[i] = right_result[i];
#    }

#    arr_hist = arr_hist.slice(0,bincount);
#    left = left.slice(0,bincount);
#    right = right.slice(0,bincount);

#    data['arr_hist'] = arr_hist;
#    data['left'] = left;
#    data['right'] = right;

#    x_range.start = oldxMin;
#    x_range.end = oldxMax;

#    src.change.emit();

#    p.change.emit();
#""")

#callback_subtractbins = CustomJS(args=dict(src=src, p=p, x_range=p.x_range), code="""

#    var oldxMin = xMin;
#    var oldxMax = xMax;

#    p.reset.emit();

#    hasWeights = $('#Weighted').data('clicked');

#    var data_res = getData();
#    var y_data = data_res[0];
#    var w_data = [];

#    if (hasWeights){
#        w_data = data_res[1];
#    }


#    bincount -= 1;

#    var hist_result = getHistBins(y_data, w_data);
#    var arr_hist_result = hist_result[0];
#    var left_result = hist_result[1];
#    var right_result = hist_result[2];

#    var data = src.data;

#    var arr_hist = [];
#    var left = [];
#    var right = [];

#    arr_hist = data['arr_hist'];
#    left = data['left'];
#    right = data['right'];

#    for (var i = 0; i < bincount; i++){
#        arr_hist[i] = arr_hist_result[i];
#        left[i] = left_result[i];
#        right[i] = right_result[i];
#    }

#    arr_hist.pop();
#    left.pop();
#    right.pop();

#    arr_hist = arr_hist.slice(0,bincount);
#    left = left.slice(0,bincount);
#    right = right.slice(0,bincount);

#    data['arr_hist'] = arr_hist;
#    data['left'] = left;
#    data['right'] = right;

#    x_range.start = oldxMin;
#    x_range.end = oldxMax;

#    src.change.emit();
#    p.change.emit();
#""")

callbackSelectHistParam = CustomJS(code="""
    let val = this.value;
    whichParam = Number(val.slice(-1));
""")

#layout
options = [];
#will label as "a0", "a1", etc

select = Select(title = "Select Model Parameter",
                options = options,
                value = "",
                )

callbackGenerateDropdown = CustomJS(args=dict(select=select), code="""
    let options = [];
    for (let m = 0; m < parameters.length - 2; m++){ //excluding slop
        options.push("a" + m.toString());
    }
    select.options = options;
""")

widgets = row(select)
layout = column(widgets, p)

p.js_on_event(MouseEnter, callbackGenerateDropdown)
p.js_on_event(MouseEnter, callback_plot)

select.js_on_event(MouseEnter, callbackGenerateDropdown)
select.js_on_change('value', callbackSelectHistParam)
select.js_on_change('value', callback_plot)

curdoc().add_root(layout) #any changes to layout will trigger on_change callbacks
curdoc().title = "Histogram of Parameter Distributions"

#show(layout)
save(layout)

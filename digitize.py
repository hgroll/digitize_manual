"""Digitize an image.

Use matplotlib to digitize an image by selecting points for a curve.
In Figure 1, four points need to be specified in the text boxes to
set the physical axis (x0, x1, y0, y1).
Activate a button (y/x origin, ymax, xmax) and pick points for
the relative axis:
- (y/x origin) corresponds to the value given in (y0, x0)
- (ymax) corresponds to (y1)
- (xmax) corresponds to (x0)
Add a temporary curve with the button (add), followed by picking points
in the plot. Points in the temporary curve can be deleted by pressing (undo).
Once the curve is completed, press (add) again. This also activates a new
temporary curve which can be repeated. Completed curves are shown in Figure 2.
The button (save) writes an output file for each curve with the following format:
    x       y
    44.949002879965235      125.34415584415586
    72.66187578112265       127.42262131174266
    100.37474868228007      132.96519589197413
    132.93737434114007      139.8934141172635
Additionally, (save) writes Figure 2 as an image in 'png' format.

Usage:
    python digitize.py -i <img_file> -o <dat_file>

@author: Herbert Groll <herbert.groll@nt.tuwien.ac.at>
@date: 2019-01-14
"""
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from matplotlib.widgets import Button, TextBox
import numpy as np
import argparse
import csv


class DigitizeScene:
    """User interface class for Digitize."""
    def __init__(self, img_file, csv_file):
        self.src_fig, self.src_ax = plt.subplots()
        self.dst_fig, self.dst_ax = plt.subplots()
        self.click_mode = None
        self.CLICK_MODE_AXIS = {'xyorigin', 'ymax', 'xmax'}
        self.buttons = dict(xyorigin=None,
                            ymax=None,
                            xmax=None,
                            add=None,
                            save=None,
                            )
        self.button_press_event = None
        self.inputs = dict(x0=None,
                           x1=None,
                           y0=None,
                           y1=None,
                           )
        self.axis_lines = dict(xyorigin=None,
                               ymax=None,
                               xmax=None)
        self.curves_lines = []
        self.temp_curve_lines = self.src_ax.plot((), '-x')[0]

        # actual data
        self.data = Digitize()

        self.img_file = img_file
        self.csv_file = csv_file
        self.scene = self.setup_scene()

    def setup_scene(self):
        """Setup scene (figure, canvas, axes, buttons, inputs)."""
        scene = dict()
        scene['img'] = mpimg.imread(self.img_file)

        scene['aximg'] = self.src_ax.imshow(scene['img'])
        scene['aximg'].axes.get_xaxis().set_visible(False)
        scene['aximg'].axes.get_yaxis().set_visible(False)
        # remove border
        for spine_key, spine in self.src_ax.spines.items():
            spine.set_visible(False)

        scene['cidpress'] = self.src_fig.canvas.mpl_connect('button_press_event', self.button_press_handler)
        scene['cidrelease'] = self.src_fig.canvas.mpl_connect('button_release_event', self.button_release_handler)

        # place input Text Boxes
        ax_y1 = self.src_fig.add_axes([0.05, 0.9, 0.05, 0.04])
        self.inputs['y1'] = TextBox(ax_y1, label='y1:', initial=str(self.data.y_phys_lim[1]))
        self.inputs['y1'].on_submit(lambda y: self.update_phys_lim(xy=1, idx=1, value_string=y))

        ax_y0 = self.src_fig.add_axes([0.05, 0.05, 0.05, 0.04])
        self.inputs['y0'] = TextBox(ax_y0, label='y0:', initial=str(self.data.y_phys_lim[0]))
        self.inputs['y0'].on_submit(lambda y: self.update_phys_lim(xy=1, idx=0, value_string=y))

        ax_x1 = self.src_fig.add_axes([0.9, 0.03, 0.05, 0.04])
        self.inputs['x1'] = TextBox(ax_x1, label='x1:', initial=str(self.data.x_phys_lim[1]))
        self.inputs['x1'].on_submit(lambda x: self.update_phys_lim(xy=0, idx=1, value_string=x))

        ax_x0 = self.src_fig.add_axes([0.15, 0.03, 0.05, 0.04])
        self.inputs['x0'] = TextBox(ax_x0, label='x0:', initial=str(self.data.x_phys_lim[0]))
        self.inputs['x0'].on_submit(lambda x: self.update_phys_lim(xy=0, idx=0, value_string=x))

        # buttons
        ax_b0 = self.src_fig.add_axes([0.15, 0.9, 0.2, 0.04])
        self.buttons['xyorigin'] = Button(ax_b0, 'y/x origin')
        self.buttons['xyorigin'].on_clicked(lambda x: self.toggle_button('xyorigin', x))

        ax_b1 = self.src_fig.add_axes([0.35, 0.9, 0.1, 0.04])
        self.buttons['ymax'] = Button(ax_b1, 'ymax')
        self.buttons['ymax'].on_clicked(lambda x: self.toggle_button('ymax', x))

        ax_b2 = self.src_fig.add_axes([0.45, 0.9, 0.1, 0.04])
        self.buttons['xmax'] = Button(ax_b2, 'xmax')
        self.buttons['xmax'].on_clicked(lambda x: self.toggle_button('xmax', x))

        ax_badd = self.src_fig.add_axes([0.65, 0.9, 0.1, 0.04])
        self.buttons['add'] = Button(ax_badd, 'add')
        self.buttons['add'].on_clicked(lambda x: self.toggle_button('add', x))

        ax_bundo = self.src_fig.add_axes([0.75, 0.9, 0.1, 0.04])
        self.buttons['undo'] = Button(ax_bundo, 'undo')
        self.buttons['undo'].on_clicked(self.undo_point)

        ax_bsave = self.src_fig.add_axes([0.85, 0.9, 0.1, 0.04])
        self.buttons['save'] = Button(ax_bsave, 'save')
        self.buttons['save'].on_clicked(self.save_traces)

        # set title
        self.set_title_dst()

        plt.show()
        return scene

    def button_press_handler(self, event):
        """Remember mouse position on button press."""
        if event.button == 1 and event.dblclick is False:
            self.button_press_event = event
        else:
            self.button_press_event = None

    def button_release_handler(self, event):
        """Do action on button release."""

        # ignore most events
        if event.inaxes is not self.src_ax:
            # ignore click not in axis
            return
        elif (event.button != 1
            or event.dblclick):
            # ignore wrong button
            return
        elif self.button_press_event is None:
            # ignore release without press
            return
        elif (self.button_press_event.x != event.x
            or self.button_press_event.y != event.y):
            # ignore rectangle movement (press, move, release)
            return
        else:
            # clear button press
            self.button_press_event = None

        # modify axis points and add curve points
        if self.click_mode in self.CLICK_MODE_AXIS:
            self.data.axis_points[self.click_mode] = (event.xdata, event.ydata)
            self.update_axis(self.click_mode)
            self.deactivate_button()
        elif self.click_mode == 'add':
            self.data.add_point((event.xdata, event.ydata))
            self.refresh_curve_line()

    def deactivate_button(self):
        """"Deselect a previously selected button."""
        old_button = self.buttons.get(self.click_mode, None)
        if old_button is not None:
            old_button.color = 'green'

        self.click_mode = None

    def activate_button(self, key):
        """Activate button."""
        new_button = self.buttons.get(key, None)
        if new_button is not None:
            new_button.color = 'red'

        self.click_mode = key

    def toggle_button(self, key, event):
        """Handle a button press."""
        if self.click_mode == 'add':
            self.data.accept_curve()
            self.prepare_next_curve_line()

        self.deactivate_button()

        if self.click_mode is not key:
            self.activate_button(key)

    def check_value(self, value):
        """Check user inputs."""
        try:
            return float(value)
        except ValueError:
            raise ValueError("Value is not convertible into floating point.")

    def update_phys_lim(self, xy=0, idx=0, value_string=0):
        """Update physical limits with user input."""
        value = self.check_value(value_string)

        if xy == 0:
            self.data.x_phys_lim[idx] = value
        else:
            self.data.y_phys_lim[idx] = value

        self.refresh_dst()

    def refresh_curve_line(self):
        """Update figure's temporary curve."""
        self.temp_curve_lines.set_data(self.data.temp_curve)
        self.src_fig.canvas.draw()

    def prepare_next_curve_line(self):
        """Accept temporary curve as curve and prepare new temporary curve."""
        self.curves_lines.append(self.temp_curve_lines)
        self.temp_curve_lines = self.src_ax.plot((), '-x')[0]
        self.refresh_dst()

    def set_title_dst(self):
        """Set title in destination figure."""
        self.dst_ax.set_title('digitized plot for {}'.format(self.img_file))

    def refresh_dst(self):
        """Refresh figure for finished curves."""
        self.dst_ax.clear()
        self.set_title_dst()
        for x, y in self.data.get_curves():
            self.dst_ax.plot(x, y)
        self.dst_fig.canvas.draw()

    def update_axis(self, key):
        """Update axis point value and visualization."""
        # if axis point already exists, remove it first
        if self.axis_lines[key] is not None:
            self.axis_lines[key].remove()

        xdata, ydata = self.data.axis_points[key]
        self.axis_lines[key] = self.src_ax.plot(xdata, ydata, 'x')[0]
        self.refresh_dst()

    def save_traces(self, event):
        """Invoke curves saving."""
        print('save traces')
        self.data.write_curves(self.csv_file)
        digi_filename = self.img_file[:-4] + "_digi.png"
        self.dst_fig.savefig(digi_filename, format='png')

    def undo_point(self, event):
        """Delete a point from temporary curve."""
        self.data.del_point()
        self.refresh_curve_line()


class Digitize:
    """Digitize class handling data."""
    def __init__(self):
        self.x_phys_lim = np.array([0, 1])
        self.y_phys_lim = np.array([0, 1])
        self.axis_points = dict(xyorigin=np.array([0, 1]),
                                ymax=np.array([0, 1]),
                                xmax=np.array([0, 1]),
                                )
        # store entered points
        self.temp_curve = ([],[])
        # store full curves
        self.curves = ([],[])

    def add_point(self, point):
        """Append a single point to temporary curve."""
        xdata, ydata = point
        self.temp_curve[0].append(xdata)
        self.temp_curve[1].append(ydata)

    def del_point(self):
        """Delete the last temporary point."""
        if (len(self.temp_curve[0]) > 0):
            self.temp_curve[0].pop(-1)
            self.temp_curve[1].pop(-1)

    def accept_curve(self):
        """Complete temporary curve and append to finished curves."""
        self.curves[0].append(np.array(self.temp_curve[0]))
        self.curves[1].append(np.array(self.temp_curve[1]))
        self.temp_curve[0].clear()
        self.temp_curve[1].clear()
        
    def get_curves(self):
        """Generator for final curves."""
        for x, y in zip(self.curves[0], self.curves[1]):
            rel_x = x - self.axis_points['xyorigin'][0]
            rel_y = y - self.axis_points['xyorigin'][1]
            scale_x = (self.x_phys_lim[1] - self.x_phys_lim[0]) \
                      / (self.axis_points['xmax'][0] - self.axis_points['xyorigin'][0])
            scale_y = (self.y_phys_lim[1] - self.y_phys_lim[0]) \
                      / (self.axis_points['ymax'][1] - self.axis_points['xyorigin'][1])

            yield (rel_x * scale_x, rel_y * scale_y)

    def write_curves(self, filename):
        """Write each curve to a file."""
        num_curves = len(self.curves[0])

        # create -0.dat, -1.dat, ... post-fix for multiple curves
        if num_curves == 1:
            self.write_csv(filename, [self.curves[0][0], self.curves[1][0]], header=['x', 'y'])
        if num_curves > 1:
            for idx, curve_xy in enumerate(zip(self.curves[0], self.curves[1])):
                postfixed_filename = "{name}-{postfix}.dat".format(name=filename[:-4], postfix=idx)
                self.write_csv(postfixed_filename,  curve_xy, header=['x', 'y'])

    @staticmethod
    def write_csv(filename, columns, header=None, my_delimiter='\t'):
        """Write multiple columns into a csv file."""
        with open(filename, 'w', newline="\n", encoding="utf-8") as csvfile:
            my_writer = csv.writer(csvfile, delimiter=my_delimiter)
            if header is not None:
                my_writer.writerow(header)

            for row in zip(*columns):
                my_writer.writerow(row)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(prog='digitize', allow_abbrev=False)
    parser.add_argument('-i', type=str, default='Nel07a.png')
    parser.add_argument('-o', type=str, default='Nel07a.dat')
    args=parser.parse_args()

    # setup scene
    digi = DigitizeScene(args.i, args.o)
"""
Copyright (C) 2018-2021 RISE Research Institute of Sweden AB

File: window.py

Author: anders.holst@ri.se

"""


# %matplotlib notebook

import matplotlib.pyplot as plt
import matplotlib as mpl

mpl.interactive(True)
mpl.rcParams['toolbar']='None'

class Wind():
    def __init__(self, name, width, height):
        self.width = width
        self.height = height
        self.fig = plt.figure(name, (width/100.0, height/100.0))
        self.trans = self.fig.transFigure.inverted()
        self.bpfid = -1
        self.brfid = -1
        self.whfid = -1
        self.mfid = -1
        self.kfid = -1
        self.buttons = {}
        self.scrolls = {}
        self.keys = {}
        self.drags = {}
        self.objs = []
        self.objtps = []
        self.bpressed = False

    def add_line(self, x1, y1, x2, y2, wdt, col):
        fr = self.trans.transform((x1, y1))
        to = self.trans.transform((x2, y2))
        obj = plt.Line2D((fr[0], to[0]), (fr[1], to[1]), linewidth=wdt*0.75, color=col)
        self.fig.add_artist(obj)
        return obj

    def add_rect(self, x1, y1, wdt, hgt, brd, fg, bg):
        fr = self.trans.transform((x1, y1))
        to = self.trans.transform((x1+wdt, y1+hgt))
        if bg:
            obj = plt.Rectangle(fr, to[0]-fr[0], to[1]-fr[1], linewidth=brd*0.75 or 0, edgecolor=fg or (0,0,0,0), facecolor=bg)
        else:
            obj = plt.Rectangle(fr, to[0]-fr[0], to[1]-fr[1], linewidth=brd*0.75 or 0, edgecolor=fg or (0,0,0,0), fill=False)
        self.fig.add_artist(obj)
        return obj

    def add_ellipse(self, x1, y1, wdt, hgt, brd, fg, bg, angle = 0):
        fr = self.trans.transform((x1, y1))
        to = self.trans.transform((x1+wdt, y1+hgt))
        if bg:
            obj = mpl.patches.Ellipse(((fr[0]+to[0])*0.5, (fr[1]+to[1])*0.5), to[0]-fr[0], to[1]-fr[1], linewidth=brd*0.75 or 0, edgecolor=fg or (0,0,0,0), facecolor=bg, angle = angle)
        else:
            obj = mpl.patches.Ellipse(fr, to[0]-fr[0], to[1]-fr[1], linewidth=brd*0.75 or 0, edgecolor=fg or (0,0,0,0), fill=False, angle = angle)
        self.fig.add_artist(obj)
        return obj

    def add_polygon(self, xylst, brd, fg, bg):
        xytr = list(map(self.trans.transform, xylst))
        if bg:
            obj = plt.Polygon(xytr, linewidth=brd*0.75 or 0, edgecolor=fg or (0,0,0,0), facecolor=bg)
        else:
            obj = plt.Polygon(xytr, linewidth=brd*0.75 or 0, edgecolor=fg or (0,0,0,0), fill=False)
        self.fig.add_artist(obj)
        return obj

    def add_text(self, x1, y1, txt, fontsz = None, align = None):
        fr = self.trans.transform((x1, y1))
        obj = plt.Text(fr[0], fr[1], txt, fontsize=fontsz or 18, ha=align or 'left')
        self.fig.add_artist(obj)
        return obj

    def locate_object(self, event, objs):
        for obj in objs:
            if obj.contains(event)[0]:
                return obj

    def locate_object_type(self, event, types):
        for obj in self.fig.artists:
            if True in types or type(obj) in types:
                if obj.contains(event)[0]:
                    return obj

    def button_press_callback(self, event):
        obj = self.locate_object(event, self.objs)
        if obj is False:
            obj = self.locate_object_type(event, self.objtps)
        func = False
        if obj is not False:
            if (event.guiEvent.num, event.guiEvent.state, obj) in self.buttons:
                func = self.buttons[(event.guiEvent.num, event.guiEvent.state, obj)]
            elif (event.guiEvent.num, -1, obj) in self.buttons:
                func = self.buttons[(event.guiEvent.num, -1, obj)]
            elif (event.guiEvent.num, event.guiEvent.state, True) in self.buttons:
                func = self.buttons[(event.guiEvent.num, event.guiEvent.state, True)]
            elif (event.guiEvent.num, -1, True) in self.buttons:
                func = self.buttons[(event.guiEvent.num, -1, True)]
        if func is not False:
            if self.bpressed is False:
                if type(func)==tuple:
                    self.bpressed = (False, func[1], func[2] if len(func)==3 else False, obj)
                    func = func[0]
                func(self, event, obj)
        else:
            if (event.guiEvent.num, event.guiEvent.state) in self.buttons:
                func = self.buttons[(event.guiEvent.num, event.guiEvent.state)]
            elif (event.guiEvent.num, -1) in self.buttons:
                func = self.buttons[(event.guiEvent.num, -1)]
            if func is not False:
                if self.bpressed is False:
                    if type(func)==tuple:
                        self.bpressed = (False, func[1], func[2] if len(func)==3 else False, False)
                        func = func[0]
                    func(self, event)
        # kolla om ev över object i objs
        # kolla om bunden: (ev, mod, obj), (ev, -1, obj), (ev, mod, T), (ev, -1, T), (ev, mod), (ev, -1)
        # om release-cb, spara undan (vad göra om redan finns? ignorera press?)
        # kör press-cb

    def button_release_callback(self, event):
        if self.bpressed is not False:
            func = self.bpressed[1]
            obj = self.bpressed[3]
            self.bpressed = False
            if obj is False:
                func(self, event)
            else:
                func(self, event, obj)

    def button_motion_callback(self, event):
        if self.bpressed is not False:
            func = self.bpressed[2]
            obj = self.bpressed[3]
            if obj is False:
                func(self, event)
            else:
                func(self, event, obj)

    def bind_button(self, button, func, obj = False):
        # func called as func(win, ev) or func(win, ev, obj) if object is not False
        if (self.bpfid == -1):
            self.bpfid = self.fig.canvas.mpl_connect('button_press_event', self.button_press_callback)
        if (self.brfid == -1):
            self.brfid = self.fig.canvas.mpl_connect('button_release_event', self.button_release_callback)
        if (self.mfid == -1 and type(func)==tuple and len(func)==3):
            self.mfid = self.fig.canvas.mpl_connect('motion_notify_event', self.button_motion_callback)
        but = button if type(button)==tuple else (button, -1)
        if obj is not False:
            but = but + (obj,)
            if obj is True or type(obj)==type:
                if not obj in self.objtps:
                    self.objtps.append(obj)
            elif not obj in self.objs:
                self.objs.append(obj)
        self.buttons[but] = func
        
    def key_press_callback(self, event):
        obj = self.locate_object(event, self.objs)
        if obj is False:
            obj = self.locate_object_type(event, self.objtps)
        func = False
        if obj is not False:
            if (event.guiEvent.keysym, event.guiEvent.state, obj) in self.keys:
                func = self.keys[(event.guiEvent.keysum, event.guiEvent.state, obj)]
            elif (event.guiEvent.keysym, -1, obj) in self.keys:
                func = self.keys[(event.guiEvent.keysym, -1, obj)]
            elif (event.guiEvent.keysym, event.guiEvent.state, True) in self.keys:
                func = self.keys[(event.guiEvent.keysym, event.guiEvent.state, True)]
            elif (event.guiEvent.keysym, -1, True) in self.keys:
                func = self.keys[(event.guiEvent.keysym, -1, True)]
        if func is not False:
            func(self, event, obj)
        else:
            if (event.guiEvent.keysym, event.guiEvent.state) in self.keys:
                func = self.keys[(event.guiEvent.keysym, event.guiEvent.state)]
            elif (event.guiEvent.keysym, -1) in self.keys:
                func = self.keys[(event.guiEvent.keysym, -1)]
            if func is not False:
                func(self, event)
        # kolla om ev över object i objs
        # kolla om bunden: (ev, mod, obj), (ev, -1, obj), (ev, mod, T), (ev, -1, T), (ev, mod), (ev, -1)
        # kör press-cb

    def bind_key(self, keysym, func, obj = False):
        # func called as func(win, ev) or func(win, ev, obj) if object is not False
        if (self.kfid == -1):
            self.kfid = self.fig.canvas.mpl_connect('key_press_event', self.key_press_callback)
        key = keysym if type(keysym)==tuple else (keysym, -1)
        if obj is not False:
            key = key + (obj,)
            if obj is True or type(obj)==type:
                if not obj in self.objtps:
                    self.objtps.append(obj)
            elif not obj in self.objs:
                self.objs.append(obj)
        self.keys[key] = func

    def scroll_callback(self, event):
        obj = self.locate_object(event, self.objs)
        if obj is False:
            obj = self.locate_object_type(event, self.objtps)
        func = False
        if obj is not False:
            if (event.button, event.guiEvent.state, obj) in self.scrolls:
                func = self.scrolls[(event.button, event.guiEvent.state, obj)]
            elif (event.button, -1, obj) in self.scrolls:
                func = self.scrolls[(event.button, -1, obj)]
            elif (event.button, event.guiEvent.state, True) in self.scrolls:
                func = self.scrolls[(event.button, event.guiEvent.state, True)]
            elif (event.button, -1, True) in self.scrolls:
                func = self.scrolls[(event.button, -1, True)]
            elif (True, -1, obj) in self.scrolls:
                func = self.scrolls[(True, -1, obj)]
            elif (True, -1, True) in self.scrolls:
                func = self.scrolls[(True, -1, True)]
        if func is not False:
            func(self, event, obj)
        else:
            if (event.button, event.guiEvent.state) in self.scrolls:
                func = self.keys[(event.button, event.guiEvent.state)]
            elif (event.button, -1) in self.keys:
                func = self.keys[(event.button, -1)]
            elif (True, -1) in self.scrolls:
                func = self.scrolls[(True, -1)]
            if func is not False:
                func(self, event)
        # kolla om ev över object i objs
        # kolla om bunden: (ev, mod, obj), (ev, -1, obj), (ev, mod, T), (ev, -1, T), (ev, mod), (ev, -1)
        # kör press-cb

    def bind_scroll(self, dir, func, obj = False):
        # dir can be "up", "down", "left", "right", or True
        # func called as func(win, ev) or func(win, ev, obj) if object is not False
        if (self.whfid == -1):
            self.whfid = self.fig.canvas.mpl_connect('scroll_event', self.scroll_callback)
        dd = dir if type(dir)==tuple else (dir, -1)
        if obj is not False:
            dir = dir + (obj,)
            if obj is True or type(obj)==type:
                if not obj in self.objtps:
                    self.objtps.append(obj)
            elif not obj in self.objs:
                self.objs.append(obj)
        self.keys[dd] = func



"""
Pmw copyright

Copyright 1997-1999 Telstra Corporation Limited,
Australia Copyright 2000-2002 Really Good Software Pty Ltd, Australia

Permission is hereby granted, free of charge, to any person obtaining
a copy of this software and associated documentation files (the
"Software"), to deal in the Software without restriction, including
without limitation the rights to use, copy, modify, merge, publish,
distribute, sublicense, and/or sell copies of the Software, and to
permit persons to whom the Software is furnished to do so, subject to
the following conditions:

The above copyright notice and this permission notice shall be
included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.


"""


import PmwColor
Color = PmwColor
del PmwColor

import PmwBlt
Blt = PmwBlt
del PmwBlt


### Loader functions:

_VERSION = '1.2'

def setversion(version):
    if version != _VERSION:
        raise ValueError('Dynamic versioning not available')

def setalphaversions(*alpha_versions):
    if alpha_versions != ():
	raise ValueError('Dynamic versioning not available')

def version(alpha = 0):
    if alpha:
        return ()
    else:
        return _VERSION

def installedversions(alpha = 0):
    if alpha:
        return ()
    else:
        return (_VERSION,)


######################################################################
### File: PmwBase.py
# Pmw megawidget base classes.

# This module provides a foundation for building megawidgets.  It
# contains the MegaArchetype class which manages component widgets and
# configuration options.  Also provided are the MegaToplevel and
# MegaWidget classes, derived from the MegaArchetype class.  The
# MegaToplevel class contains a Tkinter Toplevel widget to act as the
# container of the megawidget.  This is used as the base class of all
# megawidgets that are contained in their own top level window, such
# as a Dialog window.  The MegaWidget class contains a Tkinter Frame
# to act as the container of the megawidget.  This is used as the base
# class of all other megawidgets, such as a ComboBox or ButtonBox.
#
# Megawidgets are built by creating a class that inherits from either
# the MegaToplevel or MegaWidget class.

import os
import string
import sys
import traceback
import types
import Tkinter

# Special values used in index() methods of several megawidgets.
END = ['end']
SELECT = ['select']
DEFAULT = ['default']

# Constant used to indicate that an option can only be set by a call
# to the constructor.
INITOPT = ['initopt']
_DEFAULT_OPTION_VALUE = ['default_option_value']
_useTkOptionDb = 0

# Symbolic constants for the indexes into an optionInfo list.
_OPT_DEFAULT         = 0
_OPT_VALUE           = 1
_OPT_FUNCTION        = 2

# Stacks

_busyStack = []
    # Stack which tracks nested calls to show/hidebusycursor (called
    # either directly or from activate()/deactivate()).  Each element
    # is a dictionary containing:
    #   'newBusyWindows' :  List of windows which had busy_hold called
    #                       on them during a call to showbusycursor(). 
    #                       The corresponding call to hidebusycursor()
    #                       will call busy_release on these windows.
    #   'busyFocus' :       The blt _Busy window which showbusycursor()
    #                       set the focus to.
    #   'previousFocus' :   The focus as it was when showbusycursor()
    #                       was called.  The corresponding call to
    #                       hidebusycursor() will restore this focus if
    #                       the focus has not been changed from busyFocus.

_grabStack = []
    # Stack of grabbed windows.  It tracks calls to push/popgrab()
    # (called either directly or from activate()/deactivate()).  The
    # window on the top of the stack is the window currently with the
    # grab.  Each element is a dictionary containing:
    #   'grabWindow' :      The window grabbed by pushgrab().  The
    #                       corresponding call to popgrab() will release
    #                       the grab on this window and restore the grab
    #                       on the next window in the stack (if there is one).
    #   'globalMode' :      True if the grabWindow was grabbed with a
    #                       global grab, false if the grab was local
    #                       and 'nograb' if no grab was performed.
    #   'previousFocus' :   The focus as it was when pushgrab()
    #                       was called.  The corresponding call to
    #                       popgrab() will restore this focus.
    #   'deactivateFunction' :
    #       The function to call (usually grabWindow.deactivate) if
    #       popgrab() is called (usually from a deactivate() method)
    #       on a window which is not at the top of the stack (that is,
    #       does not have the grab or focus).  For example, if a modal
    #       dialog is deleted by the window manager or deactivated by
    #       a timer.  In this case, all dialogs above and including
    #       this one are deactivated, starting at the top of the
    #       stack.

    # Note that when dealing with focus windows, the name of the Tk
    # widget is used, since it may be the '_Busy' window, which has no
    # python instance associated with it.

#=============================================================================

# Functions used to forward methods from a class to a component.

# Fill in a flattened method resolution dictionary for a class (attributes are 
# filtered out). Flattening honours the MI method resolution rules 
# (depth-first search of bases in order). The dictionary has method names
# for keys and functions for values.
def __methodDict(cls, dict):

    # the strategy is to traverse the class in the _reverse_ of the normal
    # order, and overwrite any duplicates.
    baseList = list(cls.__bases__)
    baseList.reverse()
    
    # do bases in reverse order, so first base overrides last base
    for super in baseList:
	__methodDict(super, dict)

    # do my methods last to override base classes
    for key, value in cls.__dict__.items():
	# ignore class attributes
	if type(value) == types.FunctionType:
	    dict[key] = value

def __methods(cls):
    # Return all method names for a class.

    # Return all method names for a class (attributes are filtered
    # out).  Base classes are searched recursively.

    dict = {}
    __methodDict(cls, dict)
    return dict.keys()
	
# Function body to resolve a forwarding given the target method name and the 
# attribute name. The resulting lambda requires only self, but will forward 
# any other parameters.
__stringBody = (
    'def %(method)s(this, *args, **kw): return ' +
    'apply(this.%(attribute)s.%(method)s, args, kw)')

# Get a unique id
__counter = 0
def __unique():
  global __counter
  __counter = __counter + 1
  return str(__counter)

# Function body to resolve a forwarding given the target method name and the
# index of the resolution function. The resulting lambda requires only self, 
# but will forward any other parameters. The target instance is identified 
# by invoking the resolution function.
__funcBody = (
    'def %(method)s(this, *args, **kw): return ' +
    'apply(this.%(forwardFunc)s().%(method)s, args, kw)')

def forwardmethods(fromClass, toClass, toPart, exclude = ()):
    # Forward all methods from one class to another.

    # Forwarders will be created in fromClass to forward method
    # invocations to toClass.  The methods to be forwarded are
    # identified by flattening the interface of toClass, and excluding
    # methods identified in the exclude list.  Methods already defined
    # in fromClass, or special methods with one or more leading or
    # trailing underscores will not be forwarded.

    # For a given object of class fromClass, the corresponding toClass
    # object is identified using toPart.  This can either be a String
    # denoting an attribute of fromClass objects, or a function taking
    # a fromClass object and returning a toClass object.

    # Example:
    #     class MyClass:
    #     ...
    #         def __init__(self):
    #             ...
    #             self.__target = TargetClass()
    #             ...
    #         def findtarget(self):
    #             return self.__target
    #     forwardmethods(MyClass, TargetClass, '__target', ['dangerous1', 'dangerous2'])
    #     # ...or...
    #     forwardmethods(MyClass, TargetClass, MyClass.findtarget, 
    #             ['dangerous1', 'dangerous2'])

    # In both cases, all TargetClass methods will be forwarded from
    # MyClass except for dangerous1, dangerous2, special methods like
    # __str__, and pre-existing methods like findtarget.


    # Allow an attribute name (String) or a function to determine the instance
    if not isinstance(toPart, basestring):

	# check that it is something like a function
	if callable(toPart):

	    # If a method is passed, use the function within it
	    if hasattr(toPart, 'im_func'):
		toPart = toPart.im_func
		
	    # After this is set up, forwarders in this class will use
	    # the forwarding function. The forwarding function name is
	    # guaranteed to be unique, so that it can't be hidden by subclasses
	    forwardName = '__fwdfunc__' + __unique()
	    fromClass.__dict__[forwardName] = toPart

	# It's not a valid type
	else:
	    raise TypeError('toPart must be attribute name, function or method')

    # get the full set of candidate methods
    dict = {}
    __methodDict(toClass, dict)

    # discard special methods
    for ex in dict.keys():
	if ex[:1] == '_' or ex[-1:] == '_':
	    del dict[ex]
    # discard dangerous methods supplied by the caller
    for ex in exclude:
	if dict.has_key(ex):
	    del dict[ex]
    # discard methods already defined in fromClass
    for ex in __methods(fromClass):
	if dict.has_key(ex):
	    del dict[ex]

    for method, func in dict.items():
	d = {'method': method, 'func': func}
        if isinstance(toPart, basestring):
	    execString = \
		__stringBody % {'method' : method, 'attribute' : toPart}
	else:
	    execString = \
		__funcBody % {'forwardFunc' : forwardName, 'method' : method}

	exec execString in d

	# this creates a method
	fromClass.__dict__[method] = d[method]

#=============================================================================

def setgeometryanddeiconify(window, geom):
    # To avoid flashes on X and to position the window correctly on NT
    # (caused by Tk bugs).

    if os.name == 'nt' or \
            (os.name == 'posix' and sys.platform[:6] == 'cygwin'):
        # Require overrideredirect trick to stop window frame
        # appearing momentarily.
        redirect = window.overrideredirect()
        if not redirect:
            window.overrideredirect(1)
        window.deiconify()
        if geom is not None:
            window.geometry(geom)
        # Call update_idletasks to ensure NT moves the window to the
        # correct position it is raised.
        window.update_idletasks()
        window.tkraise()
        if not redirect:
            window.overrideredirect(0)
    else:
        if geom is not None:
            window.geometry(geom)

        # Problem!?  Which way around should the following two calls
        # go?  If deiconify() is called first then I get complaints
        # from people using the enlightenment or sawfish window
        # managers that when a dialog is activated it takes about 2
        # seconds for the contents of the window to appear.  But if
        # tkraise() is called first then I get complaints from people
        # using the twm window manager that when a dialog is activated
        # it appears in the top right corner of the screen and also
        # takes about 2 seconds to appear.

        #window.tkraise()
        # Call update_idletasks to ensure certain window managers (eg: 
        # enlightenment and sawfish) do not cause Tk to delay for
        # about two seconds before displaying window.
        #window.update_idletasks()
        #window.deiconify()

        window.deiconify()
        if window.overrideredirect():
            # The window is not under the control of the window manager
            # and so we need to raise it ourselves.
            window.tkraise()

#=============================================================================

class MegaArchetype:
    # Megawidget abstract root class.

    # This class provides methods which are inherited by classes
    # implementing useful bases (this class doesn't provide a
    # container widget inside which the megawidget can be built).

    def __init__(self, parent = None, hullClass = None):

	# Mapping from each megawidget option to a list of information
	# about the option
	#   - default value
	#   - current value
	#   - function to call when the option is initialised in the
	#     call to initialiseoptions() in the constructor or
	#     modified via configure().  If this is INITOPT, the
	#     option is an initialisation option (an option that can
	#     be set by the call to the constructor but can not be
	#     used with configure).
	# This mapping is not initialised here, but in the call to
	# defineoptions() which precedes construction of this base class.
	#
	# self._optionInfo = {}

	# Mapping from each component name to a tuple of information
	# about the component.
	#   - component widget instance
	#   - configure function of widget instance
	#   - the class of the widget (Frame, EntryField, etc)
	#   - cget function of widget instance
	#   - the name of the component group of this component, if any
	self.__componentInfo = {}

	# Mapping from alias names to the names of components or
	# sub-components.
	self.__componentAliases = {}

	# Contains information about the keywords provided to the
	# constructor.  It is a mapping from the keyword to a tuple
	# containing:
	#    - value of keyword
	#    - a boolean indicating if the keyword has been used.
	# A keyword is used if, during the construction of a megawidget,
	#    - it is defined in a call to defineoptions() or addoptions(), or
	#    - it references, by name, a component of the megawidget, or
	#    - it references, by group, at least one component
	# At the end of megawidget construction, a call is made to
	# initialiseoptions() which reports an error if there are
	# unused options given to the constructor.
        #
        # After megawidget construction, the dictionary contains
        # keywords which refer to a dynamic component group, so that
        # these components can be created after megawidget
        # construction and still use the group options given to the
        # constructor.
	#
	# self._constructorKeywords = {}

        # List of dynamic component groups.  If a group is included in
        # this list, then it not an error if a keyword argument for
        # the group is given to the constructor or to configure(), but
        # no components with this group have been created.
        # self._dynamicGroups = ()

	if hullClass is None:
	    self._hull = None
	else:
	    if parent is None:
		parent = Tkinter._default_root

	    # Create the hull.
	    self._hull = self.createcomponent('hull',
		    (), None,
		    hullClass, (parent,))
	    _hullToMegaWidget[self._hull] = self

	    if _useTkOptionDb:
		# Now that a widget has been created, query the Tk
		# option database to get the default values for the
		# options which have not been set in the call to the
		# constructor.  This assumes that defineoptions() is
		# called before the __init__().
		option_get = self.option_get
		_VALUE = _OPT_VALUE
		_DEFAULT = _OPT_DEFAULT
		for name, info in self._optionInfo.items():
		    value = info[_VALUE]
		    if value is _DEFAULT_OPTION_VALUE:
			resourceClass = string.upper(name[0]) + name[1:]
			value = option_get(name, resourceClass)
			if value != '':
			    try:
				# Convert the string to int/float/tuple, etc
				value = eval(value, {'__builtins__': {}})
			    except:
				pass
			    info[_VALUE] = value
			else:
			    info[_VALUE] = info[_DEFAULT]

    def destroy(self):
        # Clean up optionInfo in case it contains circular references
        # in the function field, such as self._settitle in class
        # MegaToplevel.

	self._optionInfo = {}
        if self._hull is not None:
            del _hullToMegaWidget[self._hull]
            self._hull.destroy()

    #======================================================================
    # Methods used (mainly) during the construction of the megawidget.

    def defineoptions(self, keywords, optionDefs, dynamicGroups = ()):
	# Create options, providing the default value and the method
	# to call when the value is changed.  If any option created by
	# base classes has the same name as one in <optionDefs>, the
	# base class's value and function will be overriden.

	# This should be called before the constructor of the base
	# class, so that default values defined in the derived class
	# override those in the base class.

	if not hasattr(self, '_constructorKeywords'):
            # First time defineoptions has been called.
	    tmp = {}
	    for option, value in keywords.items():
		tmp[option] = [value, 0]
	    self._constructorKeywords = tmp
	    self._optionInfo = {}
	    self._initialiseoptions_counter = 0
        self._initialiseoptions_counter = self._initialiseoptions_counter + 1

        if not hasattr(self, '_dynamicGroups'):
            self._dynamicGroups = ()
        self._dynamicGroups = self._dynamicGroups + tuple(dynamicGroups)
	self.addoptions(optionDefs)

    def addoptions(self, optionDefs):
	# Add additional options, providing the default value and the
	# method to call when the value is changed.  See
	# "defineoptions" for more details

	# optimisations:
	optionInfo = self._optionInfo
	optionInfo_has_key = optionInfo.has_key
	keywords = self._constructorKeywords
	keywords_has_key = keywords.has_key
	FUNCTION = _OPT_FUNCTION

	for name, default, function in optionDefs:
	    if '_' not in name:
		# The option will already exist if it has been defined
		# in a derived class.  In this case, do not override the
		# default value of the option or the callback function
		# if it is not None.
		if not optionInfo_has_key(name):
		    if keywords_has_key(name):
			value = keywords[name][0]
			optionInfo[name] = [default, value, function]
			del keywords[name]
		    else:
			if _useTkOptionDb:
			    optionInfo[name] = \
				    [default, _DEFAULT_OPTION_VALUE, function]
			else:
			    optionInfo[name] = [default, default, function]
		elif optionInfo[name][FUNCTION] is None:
		    optionInfo[name][FUNCTION] = function
	    else:
		# This option is of the form "component_option".  If this is
		# not already defined in self._constructorKeywords add it.
		# This allows a derived class to override the default value
		# of an option of a component of a base class.
		if not keywords_has_key(name):
		    keywords[name] = [default, 0]

    def createcomponent(self, componentName, componentAliases,
            componentGroup, widgetClass, *widgetArgs, **kw):
	# Create a component (during construction or later).

	if self.__componentInfo.has_key(componentName):
	    raise ValueError('Component "%s" already exists' % componentName)

	if '_' in componentName:
	    raise ValueError('Component name "%s" must not contain "_"'
                             % componentName)

	if hasattr(self, '_constructorKeywords'):
	    keywords = self._constructorKeywords
	else:
	    keywords = {}
	for alias, component in componentAliases:
	    # Create aliases to the component and its sub-components.
	    index = string.find(component, '_')
	    if index < 0:
		self.__componentAliases[alias] = (component, None)
	    else:
		mainComponent = component[:index]
		subComponent = component[(index + 1):]
		self.__componentAliases[alias] = (mainComponent, subComponent)

	    # Remove aliases from the constructor keyword arguments by
	    # replacing any keyword arguments that begin with *alias*
	    # with corresponding keys beginning with *component*.

	    alias = alias + '_'
	    aliasLen = len(alias)
	    for option in keywords.keys():
		if len(option) > aliasLen and option[:aliasLen] == alias:
		    newkey = component + '_' + option[aliasLen:]
		    keywords[newkey] = keywords[option]
		    del keywords[option]

	componentPrefix = componentName + '_'
	nameLen = len(componentPrefix)
	for option in keywords.keys():
	    if len(option) > nameLen and option[:nameLen] == componentPrefix:
		# The keyword argument refers to this component, so add
		# this to the options to use when constructing the widget.
		kw[option[nameLen:]] = keywords[option][0]
		del keywords[option]
	    else:
		# Check if this keyword argument refers to the group
		# of this component.  If so, add this to the options
		# to use when constructing the widget.  Mark the
		# keyword argument as being used, but do not remove it
		# since it may be required when creating another
		# component.
		index = string.find(option, '_')
		if index >= 0 and componentGroup == option[:index]:
		    rest = option[(index + 1):]
		    kw[rest] = keywords[option][0]
		    keywords[option][1] = 1

	if kw.has_key('pyclass'):
	    widgetClass = kw['pyclass']
	    del kw['pyclass']
	if widgetClass is None:
	    return None
        if len(widgetArgs) == 1 and isinstance(widgetArgs[0], tuple):
            # Arguments to the constructor can be specified as either
            # multiple trailing arguments to createcomponent() or as a
            # single tuple argument.
            widgetArgs = widgetArgs[0]
	widget = apply(widgetClass, widgetArgs, kw)
	componentClass = widget.__class__.__name__
	self.__componentInfo[componentName] = (widget, widget.configure,
		componentClass, widget.cget, componentGroup)

	return widget

    def destroycomponent(self, name):
	# Remove a megawidget component.

	# This command is for use by megawidget designers to destroy a
	# megawidget component.

	self.__componentInfo[name][0].destroy()
	del self.__componentInfo[name]

    def createlabel(self, parent, childCols = 1, childRows = 1):

	labelpos = self['labelpos']
	labelmargin = self['labelmargin']
	if labelpos is None:
	    return

	label = self.createcomponent('label',
		(), None,
		Tkinter.Label, (parent,))

	if labelpos[0] in 'ns':
	    # vertical layout
	    if labelpos[0] == 'n':
		row = 0
		margin = 1
	    else:
		row = childRows + 3
		margin = row - 1
	    label.grid(column=2, row=row, columnspan=childCols, sticky=labelpos)
	    parent.grid_rowconfigure(margin, minsize=labelmargin)
	else:
	    # horizontal layout
	    if labelpos[0] == 'w':
		col = 0
		margin = 1
	    else:
		col = childCols + 3
		margin = col - 1
	    label.grid(column=col, row=2, rowspan=childRows, sticky=labelpos)
	    parent.grid_columnconfigure(margin, minsize=labelmargin)

    def initialiseoptions(self, dummy = None):
        self._initialiseoptions_counter = self._initialiseoptions_counter - 1
	if self._initialiseoptions_counter == 0:
	    unusedOptions = []
	    keywords = self._constructorKeywords
	    for name in keywords.keys():
		used = keywords[name][1]
		if not used:
                    # This keyword argument has not been used.  If it
                    # does not refer to a dynamic group, mark it as
                    # unused.
                    index = string.find(name, '_')
                    if index < 0 or name[:index] not in self._dynamicGroups:
                        unusedOptions.append(name)
	    if len(unusedOptions) > 0:
		if len(unusedOptions) == 1:
		    text = 'Unknown option "'
		else:
		    text = 'Unknown options "'
		raise KeyError(text + string.join(unusedOptions, ', ')
                               + '" for ' + self.__class__.__name__)

	    # Call the configuration callback function for every option.
	    FUNCTION = _OPT_FUNCTION
	    for info in self._optionInfo.values():
		func = info[FUNCTION]
		if func is not None and func is not INITOPT:
		    func()

    #======================================================================
    # Method used to configure the megawidget.

    def configure(self, option=None, **kw):
	# Query or configure the megawidget options.
	#
	# If not empty, *kw* is a dictionary giving new
	# values for some of the options of this megawidget or its
	# components.  For options defined for this megawidget, set
	# the value of the option to the new value and call the
	# configuration callback function, if any.  For options of the
	# form <component>_<option>, where <component> is a component
	# of this megawidget, call the configure method of the
	# component giving it the new value of the option.  The
	# <component> part may be an alias or a component group name.
	#
	# If *option* is None, return all megawidget configuration
	# options and settings.  Options are returned as standard 5
	# element tuples
	#
	# If *option* is a string, return the 5 element tuple for the
	# given configuration option.

	# First, deal with the option queries.
	if len(kw) == 0:
	    # This configure call is querying the values of one or all options.
	    # Return 5-tuples:
	    #     (optionName, resourceName, resourceClass, default, value)
	    if option is None:
		rtn = {}
		for option, config in self._optionInfo.items():
		    resourceClass = string.upper(option[0]) + option[1:]
		    rtn[option] = (option, option, resourceClass,
			    config[_OPT_DEFAULT], config[_OPT_VALUE])
		return rtn
	    else:
		config = self._optionInfo[option]
		resourceClass = string.upper(option[0]) + option[1:]
		return (option, option, resourceClass, config[_OPT_DEFAULT],
			config[_OPT_VALUE])

	# optimisations:
	optionInfo = self._optionInfo
	optionInfo_has_key = optionInfo.has_key
	componentInfo = self.__componentInfo
	componentInfo_has_key = componentInfo.has_key
	componentAliases = self.__componentAliases
	componentAliases_has_key = componentAliases.has_key
	VALUE = _OPT_VALUE
	FUNCTION = _OPT_FUNCTION

	# This will contain a list of options in *kw* which
	# are known to this megawidget.
	directOptions = []

	# This will contain information about the options in
	# *kw* of the form <component>_<option>, where
	# <component> is a component of this megawidget.  It is a
	# dictionary whose keys are the configure method of each
	# component and whose values are a dictionary of options and
	# values for the component.
	indirectOptions = {}
	indirectOptions_has_key = indirectOptions.has_key

	for option, value in kw.items():
	    if optionInfo_has_key(option):
		# This is one of the options of this megawidget. 
		# Make sure it is not an initialisation option.
		if optionInfo[option][FUNCTION] is INITOPT:
		    raise KeyError('Cannot configure initialisation option "'
                                   + option + '" for '
                                   + self.__class__.__name__)
		optionInfo[option][VALUE] = value
		directOptions.append(option)
	    else:
		index = string.find(option, '_')
		if index >= 0:
		    # This option may be of the form <component>_<option>.
		    component = option[:index]
		    componentOption = option[(index + 1):]

		    # Expand component alias
		    if componentAliases_has_key(component):
			component, subComponent = componentAliases[component]
			if subComponent is not None:
			    componentOption = subComponent + '_' \
				    + componentOption

			# Expand option string to write on error
			option = component + '_' + componentOption

		    if componentInfo_has_key(component):
			# Configure the named component
			componentConfigFuncs = [componentInfo[component][1]]
		    else:
			# Check if this is a group name and configure all
			# components in the group.
			componentConfigFuncs = []
			for info in componentInfo.values():
			    if info[4] == component:
			        componentConfigFuncs.append(info[1])

                        if len(componentConfigFuncs) == 0 and \
                                component not in self._dynamicGroups:
			    raise KeyError('Unknown option "' + option + '" for '
                                           + self.__class__.__name__)

		    # Add the configure method(s) (may be more than
		    # one if this is configuring a component group)
		    # and option/value to dictionary.
		    for componentConfigFunc in componentConfigFuncs:
			if not indirectOptions_has_key(componentConfigFunc):
			    indirectOptions[componentConfigFunc] = {}
			indirectOptions[componentConfigFunc][componentOption] \
				= value
		else:
		    raise KeyError('Unknown option "' + option + '" for '
                                   + self.__class__.__name__)

	# Call the configure methods for any components.
	map(apply, indirectOptions.keys(),
		((),) * len(indirectOptions), indirectOptions.values())

	# Call the configuration callback function for each option.
	for option in directOptions:
	    info = optionInfo[option]
	    func = info[_OPT_FUNCTION]
	    if func is not None:
	      func()

    def __setitem__(self, key, value):
        apply(self.configure, (), {key: value})

    #======================================================================
    # Methods used to query the megawidget.

    def component(self, name):
	# Return a component widget of the megawidget given the
	# component's name
	# This allows the user of a megawidget to access and configure
	# widget components directly.

	# Find the main component and any subcomponents
	index = string.find(name, '_')
	if index < 0:
	    component = name
	    remainingComponents = None
	else:
	    component = name[:index]
	    remainingComponents = name[(index + 1):]

	# Expand component alias
	if self.__componentAliases.has_key(component):
	    component, subComponent = self.__componentAliases[component]
	    if subComponent is not None:
		if remainingComponents is None:
		    remainingComponents = subComponent
		else:
		    remainingComponents = subComponent + '_' \
			    + remainingComponents

	widget = self.__componentInfo[component][0]
	if remainingComponents is None:
	    return widget
	else:
	    return widget.component(remainingComponents)

    def interior(self):
	return self._hull

    def hulldestroyed(self):
	return not _hullToMegaWidget.has_key(self._hull)

    def __str__(self):
	return str(self._hull)

    def cget(self, option):
	# Get current configuration setting.

	# Return the value of an option, for example myWidget['font']. 

	if self._optionInfo.has_key(option):
	    return self._optionInfo[option][_OPT_VALUE]
	else:
	    index = string.find(option, '_')
	    if index >= 0:
		component = option[:index]
		componentOption = option[(index + 1):]

		# Expand component alias
		if self.__componentAliases.has_key(component):
		    component, subComponent = self.__componentAliases[component]
		    if subComponent is not None:
			componentOption = subComponent + '_' + componentOption

		    # Expand option string to write on error
		    option = component + '_' + componentOption

		if self.__componentInfo.has_key(component):
		    # Call cget on the component.
		    componentCget = self.__componentInfo[component][3]
		    return componentCget(componentOption)
		else:
		    # If this is a group name, call cget for one of
		    # the components in the group.
		    for info in self.__componentInfo.values():
			if info[4] == component:
			    componentCget = info[3]
			    return componentCget(componentOption)

	raise KeyError('Unknown option "' + option + '" for '
                       + self.__class__.__name__)

    __getitem__ = cget

    def isinitoption(self, option):
	return self._optionInfo[option][_OPT_FUNCTION] is INITOPT

    def options(self):
	options = []
	if hasattr(self, '_optionInfo'):
	    for option, info in self._optionInfo.items():
		isinit = info[_OPT_FUNCTION] is INITOPT
		default = info[_OPT_DEFAULT]
		options.append((option, default, isinit))
	    options.sort()
	return options

    def components(self):
	# Return a list of all components.

	# This list includes the 'hull' component and all widget subcomponents

	names = self.__componentInfo.keys()
	names.sort()
	return names

    def componentaliases(self):
	# Return a list of all component aliases.

	componentAliases = self.__componentAliases

	names = componentAliases.keys()
	names.sort()
	rtn = []
	for alias in names:
	    (mainComponent, subComponent) = componentAliases[alias]
	    if subComponent is None:
		rtn.append((alias, mainComponent))
	    else:
		rtn.append((alias, mainComponent + '_' + subComponent))
	    
	return rtn

    def componentgroup(self, name):
	return self.__componentInfo[name][4]

#=============================================================================

# The grab functions are mainly called by the activate() and
# deactivate() methods.
#
# Use pushgrab() to add a new window to the grab stack.  This
# releases the grab by the window currently on top of the stack (if
# there is one) and gives the grab and focus to the new widget.
#
# To remove the grab from the window on top of the grab stack, call
# popgrab().
#
# Use releasegrabs() to release the grab and clear the grab stack.

def pushgrab(grabWindow, globalMode, deactivateFunction):
    prevFocus = grabWindow.tk.call('focus')
    grabInfo = {
        'grabWindow' : grabWindow,
        'globalMode' : globalMode,
        'previousFocus' : prevFocus,
        'deactivateFunction' : deactivateFunction,
    }
    _grabStack.append(grabInfo)
    _grabtop()
    grabWindow.focus_set()

def popgrab(window):
    # Return the grab to the next window in the grab stack, if any.

    # If this window is not at the top of the grab stack, then it has
    # just been deleted by the window manager or deactivated by a
    # timer.  Call the deactivate method for the modal dialog above
    # this one on the stack. 
    if _grabStack[-1]['grabWindow'] != window:
        for index in range(len(_grabStack)):
            if _grabStack[index]['grabWindow'] == window:
                _grabStack[index + 1]['deactivateFunction']()
                break

    grabInfo = _grabStack[-1]
    del _grabStack[-1]

    topWidget = grabInfo['grabWindow']
    prevFocus = grabInfo['previousFocus']
    globalMode = grabInfo['globalMode']

    if globalMode != 'nograb':
        topWidget.grab_release()

    if len(_grabStack) > 0:
        _grabtop()
    if prevFocus != '':
        try:
            topWidget.tk.call('focus', prevFocus)
        except Tkinter.TclError:
            # Previous focus widget has been deleted. Set focus
            # to root window.
            Tkinter._default_root.focus_set()
    else:
        # Make sure that focus does not remain on the released widget.
        if len(_grabStack) > 0:
            topWidget = _grabStack[-1]['grabWindow']
            topWidget.focus_set()
        else:
            Tkinter._default_root.focus_set()

def grabstacktopwindow():
    if len(_grabStack) == 0:
        return None
    else:
        return _grabStack[-1]['grabWindow']

def releasegrabs():
    # Release grab and clear the grab stack.

    current = Tkinter._default_root.grab_current()
    if current is not None:
        current.grab_release()
    _grabStack[:] = []

def _grabtop():
    grabInfo = _grabStack[-1]
    topWidget = grabInfo['grabWindow']
    globalMode = grabInfo['globalMode']

    if globalMode == 'nograb':
        return

    while 1:
        try:
            if globalMode:
                topWidget.grab_set_global()
            else:
                topWidget.grab_set()
            break
        except Tkinter.TclError:
            # Another application has grab.  Keep trying until
            # grab can succeed.
            topWidget.after(100)

#=============================================================================

class MegaToplevel(MegaArchetype):

    def __init__(self, parent = None, **kw):
	# Define the options for this megawidget.
	optiondefs = (
            ('activatecommand',   None,                     None),
            ('deactivatecommand', None,                     None),
            ('master',            None,                     None),
            ('title',             None,                     self._settitle),
            ('hull_class',        self.__class__.__name__,  None),
	)
	self.defineoptions(kw, optiondefs)

	# Initialise the base class (after defining the options).
	MegaArchetype.__init__(self, parent, Tkinter.Toplevel)

	# Initialise instance.

        # Set WM_DELETE_WINDOW protocol, deleting any old callback, so
        # memory does not leak.
        if hasattr(self._hull, '_Pmw_WM_DELETE_name'):
            self._hull.tk.deletecommand(self._hull._Pmw_WM_DELETE_name)
        self._hull._Pmw_WM_DELETE_name = \
                self.register(self._userDeleteWindow, needcleanup = 0)
	self.protocol('WM_DELETE_WINDOW', self._hull._Pmw_WM_DELETE_name)

	# Initialise instance variables.

	self._firstShowing = 1
	# Used by show() to ensure window retains previous position on screen.

	# The IntVar() variable to wait on during a modal dialog.
	self._wait = None

	self._active = 0
	self._userDeleteFunc = self.destroy
	self._userModalDeleteFunc = self.deactivate

	# Check keywords and initialise options.
	self.initialiseoptions()

    def _settitle(self):
	title = self['title']
	if title is not None:
	    self.title(title)

    def userdeletefunc(self, func=None):
        if func:
	    self._userDeleteFunc = func
	else:
	    return self._userDeleteFunc

    def usermodaldeletefunc(self, func=None):
        if func:
	    self._userModalDeleteFunc = func
	else:
	    return self._userModalDeleteFunc

    def _userDeleteWindow(self):
	if self.active():
	    self._userModalDeleteFunc()
	else:
	    self._userDeleteFunc()

    def destroy(self):
	# Allow this to be called more than once.
	if _hullToMegaWidget.has_key(self._hull):
	    self.deactivate()

            # Remove circular references, so that object can get cleaned up.
            del self._userDeleteFunc
            del self._userModalDeleteFunc

            MegaArchetype.destroy(self)

    def show(self, master = None):
	if self.state() != 'normal':
	    if self._firstShowing:
		# Just let the window manager determine the window
		# position for the first time.
		geom = None
	    else:
		# Position the window at the same place it was last time.
		geom = self._sameposition()
            setgeometryanddeiconify(self, geom)

        if self._firstShowing:
            self._firstShowing = 0
        else:
            if self.transient() == '':
                self.tkraise()

        # Do this last, otherwise get flashing on NT:
        if master is not None:
            if master == 'parent':
                parent = self.winfo_parent()
                # winfo_parent() should return the parent widget, but the
                # the current version of Tkinter returns a string.
                if isinstance(parent, basestring):
                    parent = self._hull._nametowidget(parent)
                master = parent.winfo_toplevel()
            self.transient(master)

        self.focus()

    def _centreonscreen(self):
	# Centre the window on the screen.  (Actually halfway across
	# and one third down.)

        parent = self.winfo_parent()
        if isinstance(parent, basestring):
            parent = self._hull._nametowidget(parent)

        # Find size of window.
	self.update_idletasks()
        width = self.winfo_width()
        height = self.winfo_height()
        if width == 1 and height == 1:
            # If the window has not yet been displayed, its size is
            # reported as 1x1, so use requested size.
            width = self.winfo_reqwidth()
            height = self.winfo_reqheight()

        # Place in centre of screen:
	x = (self.winfo_screenwidth() - width) / 2 - parent.winfo_vrootx()
	y = (self.winfo_screenheight() - height) / 3 - parent.winfo_vrooty()
	if x < 0:
	    x = 0
	if y < 0:
	    y = 0
        return '+%d+%d' % (x, y)

    def _sameposition(self):
	# Position the window at the same place it was last time.

	geometry = self.geometry()
	index = string.find(geometry, '+')
	if index >= 0:
	    return geometry[index:]
        else:
	    return None

    def activate(self, globalMode = 0, geometry = 'centerscreenfirst'):
	if self._active:
	    raise ValueError('Window is already active')
	if self.state() == 'normal':
	    self.withdraw()

	self._active = 1

	showbusycursor()

	if self._wait is None:
	    self._wait = Tkinter.IntVar()
	self._wait.set(0)

	if geometry == 'centerscreenalways':
	    geom = self._centreonscreen()
	elif geometry == 'centerscreenfirst':
	    if self._firstShowing:
		# Centre the window the first time it is displayed.
		geom = self._centreonscreen()
	    else:
		# Position the window at the same place it was last time.
		geom = self._sameposition()
	elif geometry[:5] == 'first':
	    if self._firstShowing:
                geom = geometry[5:]
	    else:
		# Position the window at the same place it was last time.
		geom = self._sameposition()
        else:
            geom = geometry

	self._firstShowing = 0

        setgeometryanddeiconify(self, geom)

        # Do this last, otherwise get flashing on NT:
        master = self['master']
        if master is not None:
            if master == 'parent':
                parent = self.winfo_parent()
                # winfo_parent() should return the parent widget, but the
                # the current version of Tkinter returns a string.
                if isinstance(parent, basestring):
                    parent = self._hull._nametowidget(parent)
                master = parent.winfo_toplevel()
            self.transient(master)

        pushgrab(self._hull, globalMode, self.deactivate)
	command = self['activatecommand']
	if callable(command):
	    command()
	self.wait_variable(self._wait)

	return self._result

    def deactivate(self, result=None):
	if not self._active:
	    return
	self._active = 0

        # Restore the focus before withdrawing the window, since
        # otherwise the window manager may take the focus away so we
        # can't redirect it.  Also, return the grab to the next active
        # window in the stack, if any.
        popgrab(self._hull)

        command = self['deactivatecommand']
        if callable(command):
            command()

        self.withdraw()
        hidebusycursor(forceFocusRestore = 1)

        self._result = result
        self._wait.set(1)

    def active(self):
	return self._active

forwardmethods(MegaToplevel, Tkinter.Toplevel, '_hull')

#=============================================================================

class MegaWidget(MegaArchetype):
    def __init__(self, parent = None, **kw):
	# Define the options for this megawidget.
	optiondefs = (
	    ('hull_class',       self.__class__.__name__,  None),
	)
	self.defineoptions(kw, optiondefs)

	# Initialise the base class (after defining the options).
	MegaArchetype.__init__(self, parent, Tkinter.Frame)

	# Check keywords and initialise options.
	self.initialiseoptions()

forwardmethods(MegaWidget, Tkinter.Frame, '_hull')

#=============================================================================

# Public functions
#-----------------

_traceTk = 0
def tracetk(root = None, on = 1, withStackTrace = 0, file=None):
    global _withStackTrace
    global _traceTkFile
    global _traceTk

    if root is None:
        root = Tkinter._default_root

    _withStackTrace = withStackTrace
    _traceTk = on
    if on:
	if hasattr(root.tk, '__class__'):
	    # Tracing already on
	    return
	if file is None:
	    _traceTkFile = sys.stderr
	else:
	    _traceTkFile = file
	tk = _TraceTk(root.tk)
    else:
	if not hasattr(root.tk, '__class__'):
	    # Tracing already off
	    return
	tk = root.tk.getTclInterp()
    _setTkInterps(root, tk)

def showbusycursor():

    _addRootToToplevelBusyInfo()
    root = Tkinter._default_root

    busyInfo = {
        'newBusyWindows' : [],
        'previousFocus' : None,
        'busyFocus' : None,
    }
    _busyStack.append(busyInfo)

    if _disableKeyboardWhileBusy:
        # Remember the focus as it is now, before it is changed.
        busyInfo['previousFocus'] = root.tk.call('focus')

    if not _havebltbusy(root):
        # No busy command, so don't call busy hold on any windows.
        return

    for (window, winInfo) in _toplevelBusyInfo.items():
	if (window.state() != 'withdrawn' and not winInfo['isBusy']
                and not winInfo['excludeFromBusy']):
            busyInfo['newBusyWindows'].append(window)
            winInfo['isBusy'] = 1
            _busy_hold(window, winInfo['busyCursorName'])

            # Make sure that no events for the busy window get
            # through to Tkinter, otherwise it will crash in
            # _nametowidget with a 'KeyError: _Busy' if there is
            # a binding on the toplevel window.
            window.tk.call('bindtags', winInfo['busyWindow'], 'Pmw_Dummy_Tag')

            if _disableKeyboardWhileBusy:
                # Remember previous focus widget for this toplevel window
                # and set focus to the busy window, which will ignore all
                # keyboard events.
                winInfo['windowFocus'] = \
                        window.tk.call('focus', '-lastfor', window._w)
                window.tk.call('focus', winInfo['busyWindow'])
                busyInfo['busyFocus'] = winInfo['busyWindow']

    if len(busyInfo['newBusyWindows']) > 0:
        if os.name == 'nt':
            # NT needs an "update" before it will change the cursor.
            window.update()
        else:
            window.update_idletasks()

def hidebusycursor(forceFocusRestore = 0):

    # Remember the focus as it is now, before it is changed.
    root = Tkinter._default_root
    if _disableKeyboardWhileBusy:
        currentFocus = root.tk.call('focus')

    # Pop the busy info off the stack.
    busyInfo = _busyStack[-1]
    del _busyStack[-1]

    for window in busyInfo['newBusyWindows']:
        # If this window has not been deleted, release the busy cursor.
        if _toplevelBusyInfo.has_key(window):
            winInfo = _toplevelBusyInfo[window]
            winInfo['isBusy'] = 0
            _busy_release(window)

            if _disableKeyboardWhileBusy:
                # Restore previous focus window for this toplevel window,
                # but only if is still set to the busy window (it may have
                # been changed).
                windowFocusNow = window.tk.call('focus', '-lastfor', window._w)
                if windowFocusNow == winInfo['busyWindow']:
                    try:
                        window.tk.call('focus', winInfo['windowFocus'])
                    except Tkinter.TclError:
                        # Previous focus widget has been deleted. Set focus
                        # to toplevel window instead (can't leave focus on
                        # busy window).
                        window.focus_set()

    if _disableKeyboardWhileBusy:
        # Restore the focus, depending on whether the focus had changed
        # between the calls to showbusycursor and hidebusycursor.
        if forceFocusRestore or busyInfo['busyFocus'] == currentFocus:
            # The focus had not changed, so restore it to as it was before
            # the call to showbusycursor,
            previousFocus = busyInfo['previousFocus']
            if previousFocus is not None:
                try:
                    root.tk.call('focus', previousFocus)
                except Tkinter.TclError:
                    # Previous focus widget has been deleted; forget it.
                    pass
        else:
            # The focus had changed, so restore it to what it had been
            # changed to before the call to hidebusycursor.
            root.tk.call('focus', currentFocus)

def clearbusycursor():
    while len(_busyStack) > 0:
        hidebusycursor()

def setbusycursorattributes(window, **kw):
    _addRootToToplevelBusyInfo()
    for name, value in kw.items():
        if name == 'exclude':
            _toplevelBusyInfo[window]['excludeFromBusy'] = value
        elif name == 'cursorName':
            _toplevelBusyInfo[window]['busyCursorName'] = value
        else:
            raise KeyError('Unknown busycursor attribute "' + name + '"')

def _addRootToToplevelBusyInfo():
    # Include the Tk root window in the list of toplevels.  This must
    # not be called before Tkinter has had a chance to be initialised by
    # the application.

    root = Tkinter._default_root
    if root == None:
        root = Tkinter.Tk()
    if not _toplevelBusyInfo.has_key(root):
        _addToplevelBusyInfo(root)

def busycallback(command, updateFunction = None):
    if not callable(command):
	raise ValueError('cannot register non-command busy callback %s %s'
                         % (repr(command), type(command)))
    wrapper = _BusyWrapper(command, updateFunction)
    return wrapper.callback

_errorReportFile = None
_errorWindow = None

def reporterrorstofile(file = None):
    global _errorReportFile
    _errorReportFile = file

def displayerror(text):
    global _errorWindow

    if _errorReportFile is not None:
	_errorReportFile.write(text + '\n')
    else:
        # Print error on standard error as well as to error window. 
        # Useful if error window fails to be displayed, for example
        # when exception is triggered in a <Destroy> binding for root
        # window.
        sys.stderr.write(text + '\n')

	if _errorWindow is None:
	    # The error window has not yet been created.
	    _errorWindow = _ErrorWindow()

	_errorWindow.showerror(text)

_root = None
_disableKeyboardWhileBusy = 1

def initialise(
        root = None,
        size = None,
        fontScheme = None,
        useTkOptionDb = 0,
        noBltBusy = 0,
        disableKeyboardWhileBusy = None,
):
    # Remember if show/hidebusycursor should ignore keyboard events.
    global _disableKeyboardWhileBusy
    if disableKeyboardWhileBusy is not None:
        _disableKeyboardWhileBusy = disableKeyboardWhileBusy

    # Do not use blt busy command if noBltBusy is set.  Otherwise,
    # use blt busy if it is available.
    global _haveBltBusy
    if noBltBusy:
        _haveBltBusy = 0

    # Save flag specifying whether the Tk option database should be
    # queried when setting megawidget option default values.
    global _useTkOptionDb
    _useTkOptionDb = useTkOptionDb

    # If we haven't been given a root window, use the default or
    # create one.
    if root is None:
	if Tkinter._default_root is None:
	    root = Tkinter.Tk()
	else:
	    root = Tkinter._default_root

    # If this call is initialising a different Tk interpreter than the
    # last call, then re-initialise all global variables.  Assume the
    # last interpreter has been destroyed - ie:  Pmw does not (yet)
    # support multiple simultaneous interpreters.
    global _root
    if _root is not None and _root != root:
        global _busyStack
        global _errorWindow
        global _grabStack
        global _hullToMegaWidget
        global _toplevelBusyInfo
        _busyStack = []
        _errorWindow = None
        _grabStack = []
        _hullToMegaWidget = {}
        _toplevelBusyInfo = {}
    _root = root

    # Trap Tkinter Toplevel constructors so that a list of Toplevels
    # can be maintained.
    Tkinter.Toplevel.title = __TkinterToplevelTitle

    # Trap Tkinter widget destruction so that megawidgets can be
    # destroyed when their hull widget is destoyed and the list of
    # Toplevels can be pruned.
    Tkinter.Toplevel.destroy = __TkinterToplevelDestroy
    Tkinter.Widget.destroy = __TkinterWidgetDestroy

    # Modify Tkinter's CallWrapper class to improve the display of
    # errors which occur in callbacks.
    Tkinter.CallWrapper = __TkinterCallWrapper

    # Make sure we get to know when the window manager deletes the
    # root window.  Only do this if the protocol has not yet been set. 
    # This is required if there is a modal dialog displayed and the
    # window manager deletes the root window.  Otherwise the
    # application will not exit, even though there are no windows.
    if root.protocol('WM_DELETE_WINDOW') == '':
	root.protocol('WM_DELETE_WINDOW', root.destroy)

    # Set the base font size for the application and set the
    # Tk option database font resources.
    
    _font_initialise(root, size, fontScheme)

    return root

def alignlabels(widgets, sticky = None):
    if len(widgets) == 0:
    	return

    widgets[0].update_idletasks()

    # Determine the size of the maximum length label string.
    maxLabelWidth = 0
    for iwid in widgets:
	labelWidth = iwid.grid_bbox(0, 1)[2]
	if labelWidth > maxLabelWidth:
	    maxLabelWidth = labelWidth

    # Adjust the margins for the labels such that the child sites and
    # labels line up.
    for iwid in widgets:
	if sticky is not None:
	    iwid.component('label').grid(sticky=sticky)
	iwid.grid_columnconfigure(0, minsize = maxLabelWidth)
#=============================================================================

# Private routines
#-----------------
_callToTkReturned = 1
_recursionCounter = 1

class _TraceTk:
    def __init__(self, tclInterp):
        self.tclInterp = tclInterp

    def getTclInterp(self):
        return self.tclInterp

    # Calling from python into Tk.
    def call(self, *args, **kw):
        global _callToTkReturned
        global _recursionCounter

        _callToTkReturned = 0
        if len(args) == 1 and isinstance(args[0], tuple):
            argStr = str(args[0])
        else:
            argStr = str(args)
	_traceTkFile.write('CALL  TK> %d:%s%s' %
                (_recursionCounter, '  ' * _recursionCounter, argStr))
	_recursionCounter = _recursionCounter + 1
        try:
            result = apply(self.tclInterp.call, args, kw)
	except Tkinter.TclError, errorString:
            _callToTkReturned = 1
            _recursionCounter = _recursionCounter - 1
            _traceTkFile.write('\nTK ERROR> %d:%s-> %s\n' %
                    (_recursionCounter, '  ' * _recursionCounter,
                            repr(errorString)))
            if _withStackTrace:
                _traceTkFile.write('CALL  TK> stack:\n')
                traceback.print_stack()
            raise Tkinter.TclError(errorString)

        _recursionCounter = _recursionCounter - 1
        if _callToTkReturned:
            _traceTkFile.write('CALL RTN> %d:%s-> %s' %
                    (_recursionCounter, '  ' * _recursionCounter, repr(result)))
        else:
            _callToTkReturned = 1
            if result:
                _traceTkFile.write(' -> %s' % repr(result))
        _traceTkFile.write('\n')
        if _withStackTrace:
            _traceTkFile.write('CALL  TK> stack:\n')
            traceback.print_stack()

        _traceTkFile.flush()
        return result

    def __getattr__(self, key):
        return getattr(self.tclInterp, key)

def _setTkInterps(window, tk):
    window.tk = tk
    for child in window.children.values():
      _setTkInterps(child, tk)

#=============================================================================

# Functions to display a busy cursor.  Keep a list of all toplevels
# and display the busy cursor over them.  The list will contain the Tk
# root toplevel window as well as all other toplevel windows.
# Also keep a list of the widget which last had focus for each
# toplevel.

# Map from toplevel windows to
#     {'isBusy', 'windowFocus', 'busyWindow',
#         'excludeFromBusy', 'busyCursorName'}
_toplevelBusyInfo = {}

# Pmw needs to know all toplevel windows, so that it can call blt busy
# on them.  This is a hack so we get notified when a Tk topevel is
# created.  Ideally, the __init__ 'method' should be overridden, but
# it is a 'read-only special attribute'.  Luckily, title() is always
# called from the Tkinter Toplevel constructor.

def _addToplevelBusyInfo(window):
    if window._w == '.':
        busyWindow = '._Busy'
    else:
        busyWindow = window._w + '._Busy'

    _toplevelBusyInfo[window] = {
        'isBusy' : 0,
        'windowFocus' : None,
        'busyWindow' : busyWindow,
        'excludeFromBusy' : 0,
        'busyCursorName' : None,
    }

def __TkinterToplevelTitle(self, *args):
    # If this is being called from the constructor, include this
    # Toplevel in the list of toplevels and set the initial
    # WM_DELETE_WINDOW protocol to destroy() so that we get to know
    # about it.
    if not _toplevelBusyInfo.has_key(self):
        _addToplevelBusyInfo(self)
        self._Pmw_WM_DELETE_name = self.register(self.destroy, None, 0)
	self.protocol('WM_DELETE_WINDOW', self._Pmw_WM_DELETE_name)

    return apply(Tkinter.Wm.title, (self,) + args)

_haveBltBusy = None
def _havebltbusy(window):
    global _busy_hold, _busy_release, _haveBltBusy
    if _haveBltBusy is None:
        import PmwBlt
        _haveBltBusy = PmwBlt.havebltbusy(window)
        _busy_hold = PmwBlt.busy_hold
        if os.name == 'nt':
            # There is a bug in Blt 2.4i on NT where the busy window
            # does not follow changes in the children of a window.
            # Using forget works around the problem.
            _busy_release = PmwBlt.busy_forget
        else:
            _busy_release = PmwBlt.busy_release
    return _haveBltBusy

class _BusyWrapper:
    def __init__(self, command, updateFunction):
	self._command = command
	self._updateFunction = updateFunction

    def callback(self, *args):
	showbusycursor()
	rtn = apply(self._command, args)

	# Call update before hiding the busy windows to clear any
	# events that may have occurred over the busy windows.
	if callable(self._updateFunction):
	    self._updateFunction()

	hidebusycursor()
	return rtn

#=============================================================================

def drawarrow(canvas, color, direction, tag, baseOffset = 0.25, edgeOffset = 0.15):
    canvas.delete(tag)

    bw = (string.atoi(canvas['borderwidth']) + 
            string.atoi(canvas['highlightthickness']))
    width = string.atoi(canvas['width'])
    height = string.atoi(canvas['height'])

    if direction in ('up', 'down'):
        majorDimension = height
        minorDimension = width
    else:
        majorDimension = width
        minorDimension = height

    offset = round(baseOffset * majorDimension)
    if direction in ('down', 'right'):
        base = bw + offset
        apex = bw + majorDimension - offset
    else:
        base = bw + majorDimension - offset
        apex = bw + offset

    if minorDimension > 3 and minorDimension % 2 == 0:
        minorDimension = minorDimension - 1
    half = int(minorDimension * (1 - 2 * edgeOffset)) / 2
    low = round(bw + edgeOffset * minorDimension)
    middle = low + half
    high = low + 2 * half

    if direction in ('up', 'down'):
        coords = (low, base, high, base, middle, apex)
    else:
        coords = (base, low, base, high, apex, middle)
    kw = {'fill' : color, 'outline' : color, 'tag' : tag}
    apply(canvas.create_polygon, coords, kw)

#=============================================================================

# Modify the Tkinter destroy methods so that it notifies us when a Tk
# toplevel or frame is destroyed.

# A map from the 'hull' component of a megawidget to the megawidget. 
# This is used to clean up a megawidget when its hull is destroyed.
_hullToMegaWidget = {}

def __TkinterToplevelDestroy(tkWidget):
    if _hullToMegaWidget.has_key(tkWidget):
        mega = _hullToMegaWidget[tkWidget]
        try:
	    mega.destroy()
        except:
	    _reporterror(mega.destroy, ())
    else:
        # Delete the busy info structure for this toplevel (if the
        # window was created before initialise() was called, it
        # will not have any.
        if _toplevelBusyInfo.has_key(tkWidget):
            del _toplevelBusyInfo[tkWidget]
        if hasattr(tkWidget, '_Pmw_WM_DELETE_name'):
            tkWidget.tk.deletecommand(tkWidget._Pmw_WM_DELETE_name)
            del tkWidget._Pmw_WM_DELETE_name
        Tkinter.BaseWidget.destroy(tkWidget)

def __TkinterWidgetDestroy(tkWidget):
    if _hullToMegaWidget.has_key(tkWidget):
        mega = _hullToMegaWidget[tkWidget]
        try:
	    mega.destroy()
        except:
	    _reporterror(mega.destroy, ())
    else:
        Tkinter.BaseWidget.destroy(tkWidget)

#=============================================================================

# Add code to Tkinter to improve the display of errors which occur in
# callbacks.

class __TkinterCallWrapper:
    def __init__(self, func, subst, widget):
	self.func = func
	self.subst = subst
	self.widget = widget

    # Calling back from Tk into python.
    def __call__(self, *args):
	try:
	    if self.subst:
		args = apply(self.subst, args)
            if _traceTk:
                if not _callToTkReturned:
                    _traceTkFile.write('\n')
                if hasattr(self.func, 'im_class'):
                    name = self.func.im_class.__name__ + '.' + \
                        self.func.__name__
                else:
                    name = self.func.__name__
                if len(args) == 1 and hasattr(args[0], 'type'):
                    # The argument to the callback is an event.
                    eventName = _eventTypeToName[string.atoi(args[0].type)]
                    if eventName in ('KeyPress', 'KeyRelease',):
                        argStr = '(%s %s Event: %s)' % \
                            (eventName, args[0].keysym, args[0].widget)
                    else:
                        argStr = '(%s Event, %s)' % (eventName, args[0].widget)
                else:
                    argStr = str(args)
                _traceTkFile.write('CALLBACK> %d:%s%s%s\n' %
                    (_recursionCounter, '  ' * _recursionCounter, name, argStr))
                _traceTkFile.flush()
	    return apply(self.func, args)
	except SystemExit, msg:
	    raise SystemExit(msg)
	except:
	    _reporterror(self.func, args)

_eventTypeToName = {
    2 : 'KeyPress',         15 : 'VisibilityNotify',   28 : 'PropertyNotify',
    3 : 'KeyRelease',       16 : 'CreateNotify',       29 : 'SelectionClear',
    4 : 'ButtonPress',      17 : 'DestroyNotify',      30 : 'SelectionRequest',
    5 : 'ButtonRelease',    18 : 'UnmapNotify',        31 : 'SelectionNotify',
    6 : 'MotionNotify',     19 : 'MapNotify',          32 : 'ColormapNotify',
    7 : 'EnterNotify',      20 : 'MapRequest',         33 : 'ClientMessage',
    8 : 'LeaveNotify',      21 : 'ReparentNotify',     34 : 'MappingNotify',
    9 : 'FocusIn',          22 : 'ConfigureNotify',    35 : 'VirtualEvents',
    10 : 'FocusOut',        23 : 'ConfigureRequest',   36 : 'ActivateNotify',
    11 : 'KeymapNotify',    24 : 'GravityNotify',      37 : 'DeactivateNotify',
    12 : 'Expose',          25 : 'ResizeRequest',      38 : 'MouseWheelEvent',
    13 : 'GraphicsExpose',  26 : 'CirculateNotify',
    14 : 'NoExpose',        27 : 'CirculateRequest',
}

def _reporterror(func, args):
    # Fetch current exception values.
    exc_type, exc_value, exc_traceback = sys.exc_info()

    # Give basic information about the callback exception.
    if type(exc_type) == types.ClassType:
	# Handle python 1.5 class exceptions.
	exc_type = exc_type.__name__
    msg = exc_type + ' Exception in Tk callback\n'
    msg = msg + '  Function: %s (type: %s)\n' % (repr(func), type(func))
    msg = msg + '  Args: %s\n' % str(args)

    if isinstance(args, tuple) and len(args) > 0 and \
	    hasattr(args[0], 'type'):
        eventArg = 1
    else:
        eventArg = 0

    # If the argument to the callback is an event, add the event type.
    if eventArg:
	eventNum = string.atoi(args[0].type)
        if eventNum in _eventTypeToName.keys():
            msg = msg + '  Event type: %s (type num: %d)\n' % \
                    (_eventTypeToName[eventNum], eventNum)
        else:
            msg = msg + '  Unknown event type (type num: %d)\n' % eventNum

    # Add the traceback.
    msg = msg + 'Traceback (innermost last):\n'
    for tr in traceback.extract_tb(exc_traceback):
	msg = msg + '  File "%s", line %s, in %s\n' % (tr[0], tr[1], tr[2])
	msg = msg + '    %s\n' % tr[3]
    msg = msg + '%s: %s\n' % (exc_type, exc_value)

    # If the argument to the callback is an event, add the event contents.
    if eventArg:
	msg = msg + '\n================================================\n'
	msg = msg + '  Event contents:\n'
	keys = args[0].__dict__.keys()
	keys.sort()
	for key in keys:
	    msg = msg + '    %s: %s\n' % (key, args[0].__dict__[key])

    clearbusycursor()
    try:
	displayerror(msg)
    except:
        pass

class _ErrorWindow:
    def __init__(self):

	self._errorQueue = []
	self._errorCount = 0
	self._open = 0
        self._firstShowing = 1

	# Create the toplevel window
	self._top = Tkinter.Toplevel()
	self._top.protocol('WM_DELETE_WINDOW', self._hide)
	self._top.title('Error in background function')
	self._top.iconname('Background error')

	# Create the text widget and scrollbar in a frame
	upperframe = Tkinter.Frame(self._top)

	scrollbar = Tkinter.Scrollbar(upperframe, orient='vertical')
	scrollbar.pack(side = 'right', fill = 'y')

	self._text = Tkinter.Text(upperframe, yscrollcommand=scrollbar.set)
	self._text.pack(fill = 'both', expand = 1)
	scrollbar.configure(command=self._text.yview)

	# Create the buttons and label in a frame
	lowerframe = Tkinter.Frame(self._top)

	ignore = Tkinter.Button(lowerframe,
	        text = 'Ignore remaining errors', command = self._hide)
	ignore.pack(side='left')

	self._nextError = Tkinter.Button(lowerframe,
	        text = 'Show next error', command = self._next)
	self._nextError.pack(side='left')

	self._label = Tkinter.Label(lowerframe, relief='ridge')
	self._label.pack(side='left', fill='x', expand=1)

	# Pack the lower frame first so that it does not disappear
	# when the window is resized.
	lowerframe.pack(side = 'bottom', fill = 'x')
	upperframe.pack(side = 'bottom', fill = 'both', expand = 1)

    def showerror(self, text):
	if self._open:
	    self._errorQueue.append(text)
	else:
	    self._display(text)
	    self._open = 1

	# Display the error window in the same place it was before.
	if self._top.state() == 'normal':
	    # If update_idletasks is not called here, the window may
	    # be placed partially off the screen.  Also, if it is not
	    # called and many errors are generated quickly in
	    # succession, the error window may not display errors
	    # until the last one is generated and the interpreter
	    # becomes idle.
	    # XXX: remove this, since it causes omppython to go into an
	    # infinite loop if an error occurs in an omp callback.
	    # self._top.update_idletasks()

            pass
	else:
	    if self._firstShowing:
		geom = None
	    else:
                geometry = self._top.geometry()
                index = string.find(geometry, '+')
                if index >= 0:
                    geom = geometry[index:]
                else:
                    geom = None
            setgeometryanddeiconify(self._top, geom)

        if self._firstShowing:
            self._firstShowing = 0
        else:
            self._top.tkraise()

        self._top.focus()

	self._updateButtons()

	# Release any grab, so that buttons in the error window work.
        releasegrabs()

    def _hide(self):
	self._errorCount = self._errorCount + len(self._errorQueue)
	self._errorQueue = []
	self._top.withdraw()
	self._open = 0

    def _next(self):
	# Display the next error in the queue. 

	text = self._errorQueue[0]
	del self._errorQueue[0]

	self._display(text)
	self._updateButtons()

    def _display(self, text):
	self._errorCount = self._errorCount + 1
	text = 'Error: %d\n%s' % (self._errorCount, text)
	self._text.delete('1.0', 'end')
	self._text.insert('end', text)

    def _updateButtons(self):
	numQueued = len(self._errorQueue)
	if numQueued > 0:
	    self._label.configure(text='%d more errors' % numQueued)
	    self._nextError.configure(state='normal')
	else:
	    self._label.configure(text='No more errors')
	    self._nextError.configure(state='disabled')

######################################################################
### File: PmwDialog.py
# Based on iwidgets2.2.0/dialog.itk and iwidgets2.2.0/dialogshell.itk code.

# Convention:
#   Each dialog window should have one of these as the rightmost button:
#     Close         Close a window which only displays information.
#     Cancel        Close a window which may be used to change the state of
#                   the application.

import sys
import types
import Tkinter


# A Toplevel with a ButtonBox and child site.

class Dialog(MegaToplevel):
    def __init__(self, parent = None, **kw):

	# Define the megawidget options.
	
	optiondefs = (
	    ('buttonbox_hull_borderwidth',   1,         None),
	    ('buttonbox_hull_relief',        'raised',  None),
	    ('buttonboxpos',                 's',       INITOPT),
	    ('buttons',                      ('OK',),   self._buttons),
	    ('command',                      None,      None),
	    ('dialogchildsite_borderwidth',  1,         None),
	    ('dialogchildsite_relief',       'raised',  None),
	    ('defaultbutton',                None,      self._defaultButton),
            ('master',                       'parent',  None),
	    ('separatorwidth',               0,         INITOPT),
	)
	self.defineoptions(kw, optiondefs)

	# Initialise the base class (after defining the options).
	MegaToplevel.__init__(self, parent)

	# Create the components.

	oldInterior = MegaToplevel.interior(self)

	# Set up pack options according to the position of the button box.
        pos = self['buttonboxpos']
	if pos not in 'nsew':
	    raise ValueError('bad buttonboxpos option "%s":  should be n, s, e, or w'
                             % pos)

	if pos in 'ns':
	    orient = 'horizontal'
	    fill = 'x'
	    if pos == 'n':
	        side = 'top'
	    else:
	        side = 'bottom'
	else:
	    orient = 'vertical'
	    fill = 'y'
	    if pos == 'w':
	        side = 'left'
	    else:
	        side = 'right'

	# Create the button box.
	self._buttonBox = self.createcomponent('buttonbox',
		(), None,
		ButtonBox, (oldInterior,), orient = orient)
	self._buttonBox.pack(side = side, fill = fill)

	# Create the separating line.
	width = self['separatorwidth']
	if width > 0:
	    self._separator = self.createcomponent('separator',
		    (), None,
		    Tkinter.Frame, (oldInterior,), relief = 'sunken',
		    height = width, width = width, borderwidth = width / 2)
	    self._separator.pack(side = side, fill = fill)
	
	# Create the child site.
	self.__dialogChildSite = self.createcomponent('dialogchildsite',
		(), None,
		Tkinter.Frame, (oldInterior,))
	self.__dialogChildSite.pack(side=side, fill='both', expand=1)

	self.oldButtons = ()
	self.oldDefault = None

	self.bind('<Return>', self._invokeDefault)
        self.userdeletefunc(self._doCommand)
        self.usermodaldeletefunc(self._doCommand)
	
	# Check keywords and initialise options.
	self.initialiseoptions()

    def interior(self):
	return self.__dialogChildSite

    def invoke(self, index = DEFAULT):
	return self._buttonBox.invoke(index)

    def _invokeDefault(self, event):
	try:
	    self._buttonBox.index(DEFAULT)
	except ValueError:
	    return
	self._buttonBox.invoke()

    def _doCommand(self, name = None):
        if name is not None and self.active() and \
                grabstacktopwindow() != self.component('hull'):
            # This is a modal dialog but is not on the top of the grab
            # stack (ie:  should not have the grab), so ignore this
            # event.  This seems to be a bug in Tk and may occur in
            # nested modal dialogs.
            #
            # An example is the PromptDialog demonstration.  To
            # trigger the problem, start the demo, then move the mouse
            # to the main window, hit <TAB> and then <TAB> again.  The
            # highlight border of the "Show prompt dialog" button
            # should now be displayed.  Now hit <SPACE>, <RETURN>,
            # <RETURN> rapidly several times.  Eventually, hitting the
            # return key invokes the password dialog "OK" button even
            # though the confirm dialog is active (and therefore
            # should have the keyboard focus).  Observed under Solaris
            # 2.5.1, python 1.5.2 and Tk8.0.

            # TODO:  Give focus to the window on top of the grabstack.
            return

	command = self['command']
	if callable(command):
	    return command(name)
	else:
	    if self.active():
	        self.deactivate(name)
	    else:
	        self.withdraw()

    def _buttons(self):
	buttons = self['buttons']
	if not isinstance(buttons, (tuple, list)):
	    raise ValueError('bad buttons option "%s": should be a tuple'
                             % str(buttons))
	if self.oldButtons == buttons:
	  return

	self.oldButtons = buttons

	for index in range(self._buttonBox.numbuttons()):
	    self._buttonBox.delete(0)
	for name in buttons:
	    self._buttonBox.add(name,
		command=lambda self=self, name=name: self._doCommand(name))

	if len(buttons) > 0:
	    defaultbutton = self['defaultbutton']
	    if defaultbutton is None:
		self._buttonBox.setdefault(None)
	    else:
		try:
		    self._buttonBox.index(defaultbutton)
		except ValueError:
		    pass
		else:
		    self._buttonBox.setdefault(defaultbutton)
	self._buttonBox.alignbuttons()

    def _defaultButton(self):
	defaultbutton = self['defaultbutton']
	if self.oldDefault == defaultbutton:
	  return

	self.oldDefault = defaultbutton

	if len(self['buttons']) > 0:
	    if defaultbutton is None:
		self._buttonBox.setdefault(None)
	    else:
		try:
		    self._buttonBox.index(defaultbutton)
		except ValueError:
		    pass
		else:
		    self._buttonBox.setdefault(defaultbutton)

######################################################################
### File: PmwTimeFuncs.py
# Functions for dealing with dates and times.

import re
import string

def timestringtoseconds(text, separator = ':'):
  inputList = string.split(string.strip(text), separator)
  if len(inputList) != 3:
    raise ValueError('invalid value: ' + text)

  sign = 1
  if len(inputList[0]) > 0 and inputList[0][0] in ('+', '-'):
    if inputList[0][0] == '-':
      sign = -1
    inputList[0] = inputList[0][1:]

  if re.search('[^0-9]', string.join(inputList, '')) is not None:
    raise ValueError('invalid value: ' + text)

  hour = string.atoi(inputList[0])
  minute = string.atoi(inputList[1])
  second = string.atoi(inputList[2])

  if minute >= 60 or second >= 60:
    raise ValueError('invalid value: ' + text)
  return sign * (hour * 60 * 60 + minute * 60 + second)

_year_pivot = 50
_century = 2000

def setyearpivot(pivot, century = None):
    global _year_pivot
    global _century
    oldvalues = (_year_pivot, _century)
    _year_pivot = pivot
    if century is not None:
	_century = century
    return oldvalues

def datestringtojdn(text, format = 'ymd', separator = '/'):
  inputList = string.split(string.strip(text), separator)
  if len(inputList) != 3:
    raise ValueError('invalid value: ' + text)

  if re.search('[^0-9]', string.join(inputList, '')) is not None:
    raise ValueError('invalid value: ' + text)
  formatList = list(format)
  day = string.atoi(inputList[formatList.index('d')])
  month = string.atoi(inputList[formatList.index('m')])
  year = string.atoi(inputList[formatList.index('y')])

  if _year_pivot is not None:
    if year >= 0 and year < 100:
      if year <= _year_pivot:
	year = year + _century
      else:
	year = year + _century - 100

  jdn = ymdtojdn(year, month, day)
  if jdntoymd(jdn) != (year, month, day):
    raise ValueError('invalid value: ' + text)
  return jdn

def _cdiv(a, b):
    # Return a / b as calculated by most C language implementations,
    # assuming both a and b are integers.

    if a * b > 0:
	return a / b
    else:
	return -(abs(a) / abs(b))

def ymdtojdn(year, month, day, julian = -1, papal = 1):

    # set Julian flag if auto set
    if julian < 0:
	if papal:                          # Pope Gregory XIII's decree
	    lastJulianDate = 15821004L     # last day to use Julian calendar
	else:                              # British-American usage
	    lastJulianDate = 17520902L     # last day to use Julian calendar

	julian = ((year * 100L) + month) * 100 + day  <=  lastJulianDate

    if year < 0:
	# Adjust BC year
	year = year + 1

    if julian:
	return 367L * year - _cdiv(7 * (year + 5001L + _cdiv((month - 9), 7)), 4) + \
	    _cdiv(275 * month, 9) + day + 1729777L
    else:
	return (day - 32076L) + \
	    _cdiv(1461L * (year + 4800L + _cdiv((month - 14), 12)), 4) + \
	    _cdiv(367 * (month - 2 - _cdiv((month - 14), 12) * 12), 12) - \
	    _cdiv((3 * _cdiv((year + 4900L + _cdiv((month - 14), 12)), 100)), 4) + \
	    1            # correction by rdg

def jdntoymd(jdn, julian = -1, papal = 1):

    # set Julian flag if auto set
    if julian < 0:
	if papal:                          # Pope Gregory XIII's decree
	    lastJulianJdn = 2299160L       # last jdn to use Julian calendar
	else:                              # British-American usage
	    lastJulianJdn = 2361221L       # last jdn to use Julian calendar

	julian = (jdn <= lastJulianJdn);

    x = jdn + 68569L
    if julian:
	x = x + 38
	daysPer400Years = 146100L
	fudgedDaysPer4000Years = 1461000L + 1
    else:
	daysPer400Years = 146097L
	fudgedDaysPer4000Years = 1460970L + 31

    z = _cdiv(4 * x, daysPer400Years)
    x = x - _cdiv((daysPer400Years * z + 3), 4)
    y = _cdiv(4000 * (x + 1), fudgedDaysPer4000Years)
    x = x - _cdiv(1461 * y, 4) + 31
    m = _cdiv(80 * x, 2447)
    d = x - _cdiv(2447 * m, 80)
    x = _cdiv(m, 11)
    m = m + 2 - 12 * x
    y = 100 * (z - 49) + y + x

    # Convert from longs to integers.
    yy = int(y)
    mm = int(m)
    dd = int(d)

    if yy <= 0:
	# Adjust BC years.
	    yy = yy - 1

    return (yy, mm, dd)

def stringtoreal(text, separator = '.'):
    if separator != '.':
	if string.find(text, '.') >= 0:
	    raise ValueError('invalid value: ' + text)
	index = string.find(text, separator)
	if index >= 0:
	    text = text[:index] + '.' + text[index + 1:]
    return string.atof(text)

######################################################################
### File: PmwBalloon.py
import os
import string
import Tkinter


class Balloon(MegaToplevel):
    def __init__(self, parent = None, **kw):

	# Define the megawidget options.
	optiondefs = (
            ('initwait',                 500,            None), # milliseconds
            ('label_background',         'lightyellow',  None),
            ('label_foreground',         'black',        None),
            ('label_justify',            'left',         None),
            ('master',                   'parent',       None),
            ('relmouse',                 'none',         self._relmouse),
            ('state',                    'both',         self._state),
            ('statuscommand',            None,           None),
            ('xoffset',                  20,             None), # pixels
            ('yoffset',                  1,              None), # pixels
	    ('hull_highlightthickness',  1,              None),
	    ('hull_highlightbackground', 'black',        None),
	)
	self.defineoptions(kw, optiondefs)

	# Initialise the base class (after defining the options).
	MegaToplevel.__init__(self, parent)

	self.withdraw()
	self.overrideredirect(1)

	# Create the components.
	interior = self.interior()
	self._label = self.createcomponent('label',
		(), None,
		Tkinter.Label, (interior,))
	self._label.pack()

        # The default hull configuration options give a black border
        # around the balloon, but avoids a black 'flash' when the
        # balloon is deiconified, before the text appears.
        if not kw.has_key('hull_background'):
            self.configure(hull_background = \
                    str(self._label.cget('background')))

	# Initialise instance variables.
	self._timer = None

        # The widget or item that is currently triggering the balloon. 
        # It is None if the balloon is not being displayed.  It is a
        # one-tuple if the balloon is being displayed in response to a
        # widget binding (value is the widget).  It is a two-tuple if
        # the balloon is being displayed in response to a canvas or
        # text item binding (value is the widget and the item).
        self._currentTrigger = None
	
	# Check keywords and initialise options.
	self.initialiseoptions()

    def destroy(self):
	if self._timer is not None:
	    self.after_cancel(self._timer)
	    self._timer = None
	MegaToplevel.destroy(self)

    def bind(self, widget, balloonHelp, statusHelp = None):

        # If a previous bind for this widget exists, remove it.
        self.unbind(widget)

	if balloonHelp is None and statusHelp is None:
	    return

	if statusHelp is None:
	    statusHelp = balloonHelp
	enterId = widget.bind('<Enter>', 
		lambda event, self = self, w = widget,
			sHelp = statusHelp, bHelp = balloonHelp:
				self._enter(event, w, sHelp, bHelp, 0))

        # Set Motion binding so that if the pointer remains at rest
        # within the widget until the status line removes the help and
        # then the pointer moves again, then redisplay the help in the
        # status line.
	# Note:  The Motion binding only works for basic widgets, and
        # the hull of megawidgets but not for other megawidget components.
	motionId = widget.bind('<Motion>', 
		lambda event = None, self = self, statusHelp = statusHelp:
			self.showstatus(statusHelp))

	leaveId = widget.bind('<Leave>', self._leave)
	buttonId = widget.bind('<ButtonPress>', self._buttonpress)

        # Set Destroy binding so that the balloon can be withdrawn and
        # the timer can be cancelled if the widget is destroyed.
	destroyId = widget.bind('<Destroy>', self._destroy)

        # Use the None item in the widget's private Pmw dictionary to
        # store the widget's bind callbacks, for later clean up.
        if not hasattr(widget, '_Pmw_BalloonBindIds'):
            widget._Pmw_BalloonBindIds = {}
        widget._Pmw_BalloonBindIds[None] = \
                (enterId, motionId, leaveId, buttonId, destroyId)

    def unbind(self, widget):
        if hasattr(widget, '_Pmw_BalloonBindIds'):
            if widget._Pmw_BalloonBindIds.has_key(None):
                (enterId, motionId, leaveId, buttonId, destroyId) = \
                        widget._Pmw_BalloonBindIds[None]
                # Need to pass in old bindings, so that Tkinter can
                # delete the commands.  Otherwise, memory is leaked.
                widget.unbind('<Enter>', enterId)
                widget.unbind('<Motion>', motionId)
                widget.unbind('<Leave>', leaveId)
                widget.unbind('<ButtonPress>', buttonId)
                widget.unbind('<Destroy>', destroyId)
                del widget._Pmw_BalloonBindIds[None]

        if self._currentTrigger is not None and len(self._currentTrigger) == 1:
            # The balloon is currently being displayed and the current
            # trigger is a widget.
            triggerWidget = self._currentTrigger[0]
            if triggerWidget == widget:
                if self._timer is not None:
                    self.after_cancel(self._timer)
                    self._timer = None
                self.withdraw()
                self.clearstatus()
                self._currentTrigger = None

    def tagbind(self, widget, tagOrItem, balloonHelp, statusHelp = None):

        # If a previous bind for this widget's tagOrItem exists, remove it.
        self.tagunbind(widget, tagOrItem)

	if balloonHelp is None and statusHelp is None:
	    return

	if statusHelp is None:
	    statusHelp = balloonHelp
	enterId = widget.tag_bind(tagOrItem, '<Enter>', 
		lambda event, self = self, w = widget,
			sHelp = statusHelp, bHelp = balloonHelp:
				self._enter(event, w, sHelp, bHelp, 1))
	motionId = widget.tag_bind(tagOrItem, '<Motion>', 
		lambda event = None, self = self, statusHelp = statusHelp:
			self.showstatus(statusHelp))
	leaveId = widget.tag_bind(tagOrItem, '<Leave>', self._leave)
	buttonId = widget.tag_bind(tagOrItem, '<ButtonPress>', self._buttonpress)

        # Use the tagOrItem item in the widget's private Pmw dictionary to
        # store the tagOrItem's bind callbacks, for later clean up.
        if not hasattr(widget, '_Pmw_BalloonBindIds'):
            widget._Pmw_BalloonBindIds = {}
        widget._Pmw_BalloonBindIds[tagOrItem] = \
                (enterId, motionId, leaveId, buttonId)

    def tagunbind(self, widget, tagOrItem):
        if hasattr(widget, '_Pmw_BalloonBindIds'):
            if widget._Pmw_BalloonBindIds.has_key(tagOrItem):
                (enterId, motionId, leaveId, buttonId) = \
                        widget._Pmw_BalloonBindIds[tagOrItem]
                widget.tag_unbind(tagOrItem, '<Enter>', enterId)
                widget.tag_unbind(tagOrItem, '<Motion>', motionId)
                widget.tag_unbind(tagOrItem, '<Leave>', leaveId)
                widget.tag_unbind(tagOrItem, '<ButtonPress>', buttonId)
                del widget._Pmw_BalloonBindIds[tagOrItem]

        if self._currentTrigger is None:
            # The balloon is not currently being displayed.
            return

        if len(self._currentTrigger) == 1:
            # The current trigger is a widget.
            return

        if len(self._currentTrigger) == 2:
            # The current trigger is a canvas item.
            (triggerWidget, triggerItem) = self._currentTrigger
            if triggerWidget == widget and triggerItem == tagOrItem:
                if self._timer is not None:
                    self.after_cancel(self._timer)
                    self._timer = None
                self.withdraw()
                self.clearstatus()
                self._currentTrigger = None
        else: # The current trigger is a text item.
            (triggerWidget, x, y) = self._currentTrigger
            if triggerWidget == widget:
                currentPos = widget.index('@%d,%d' % (x, y))
                currentTags = widget.tag_names(currentPos)
                if tagOrItem in currentTags:
                    if self._timer is not None:
                        self.after_cancel(self._timer)
                        self._timer = None
                    self.withdraw()
                    self.clearstatus()
                    self._currentTrigger = None

    def showstatus(self, statusHelp):
	if self['state'] in ('status', 'both'):
	    cmd = self['statuscommand']
	    if callable(cmd):
		cmd(statusHelp)

    def clearstatus(self):
        self.showstatus(None)

    def _state(self):
	if self['state'] not in ('both', 'balloon', 'status', 'none'):
	    raise ValueError('bad state option ' + repr(self['state'])
                             + ": should be one of 'both', 'balloon', "
                             + "'status' or 'none'")

    def _relmouse(self):
	if self['relmouse'] not in ('both', 'x', 'y', 'none'):
	    raise ValueError("bad relmouse option %s: should be one of 'both', 'x', "
                             "'y' or 'none'" % repr(self['relmouse']))

    def _enter(self, event, widget, statusHelp, balloonHelp, isItem):

        # Do not display balloon if mouse button is pressed.  This
        # will only occur if the button was pressed inside a widget,
        # then the mouse moved out of and then back into the widget,
        # with the button still held down.  The number 0x1f00 is the
        # button mask for the 5 possible buttons in X.
        buttonPressed = (event.state & 0x1f00) != 0

	if not buttonPressed and balloonHelp is not None and \
                self['state'] in ('balloon', 'both'):
	    if self._timer is not None:
		self.after_cancel(self._timer)
		self._timer = None

	    self._timer = self.after(self['initwait'], 
		    lambda self = self, widget = widget, help = balloonHelp,
			    isItem = isItem:
			    self._showBalloon(widget, help, isItem))

        if isItem:
            if hasattr(widget, 'canvasx'):
		# The widget is a canvas.
                item = widget.find_withtag('current')
                if len(item) > 0:
                    item = item[0]
                else:
                    item = None
                self._currentTrigger = (widget, item)
            else:
		# The widget is a text widget.
                self._currentTrigger = (widget, event.x, event.y)
        else:
            self._currentTrigger = (widget,)

	self.showstatus(statusHelp)

    def _leave(self, event):
	if self._timer is not None:
	    self.after_cancel(self._timer)
	    self._timer = None
	self.withdraw()
	self.clearstatus()
        self._currentTrigger = None

    def _destroy(self, event):

        # Only withdraw the balloon and cancel the timer if the widget
        # being destroyed is the widget that triggered the balloon. 
        # Note that in a Tkinter Destroy event, the widget field is a
        # string and not a widget as usual.

        if self._currentTrigger is None:
            # The balloon is not currently being displayed
            return

        if len(self._currentTrigger) == 1:
            # The current trigger is a widget (not an item)
            triggerWidget = self._currentTrigger[0]
            if str(triggerWidget) == event.widget:
                if self._timer is not None:
                    self.after_cancel(self._timer)
                    self._timer = None
                self.withdraw()
                self.clearstatus()
                self._currentTrigger = None

    def _buttonpress(self, event):
	if self._timer is not None:
	    self.after_cancel(self._timer)
	    self._timer = None
	self.withdraw()
        self._currentTrigger = None

    def _showBalloon(self, widget, balloonHelp, isItem):

	self._label.configure(text = balloonHelp)

        # First, display the balloon offscreen to get dimensions.
        screenWidth = self.winfo_screenwidth()
        screenHeight = self.winfo_screenheight()
        self.geometry('+%d+0' % (screenWidth + 1))
        self.update_idletasks()

	if isItem:
            # Get the bounding box of the current item.
            bbox = widget.bbox('current')
            if bbox is None:
                # The item that triggered the balloon has disappeared,
                # perhaps by a user's timer event that occured between
                # the <Enter> event and the 'initwait' timer calling
                # this method.
                return

	    # The widget is either a text or canvas.  The meaning of
	    # the values returned by the bbox method is different for
	    # each, so use the existence of the 'canvasx' method to
	    # distinguish between them.
	    if hasattr(widget, 'canvasx'):
		# The widget is a canvas.  Place balloon under canvas
                # item.  The positions returned by bbox are relative
                # to the entire canvas, not just the visible part, so
                # need to convert to window coordinates.
                leftrel = bbox[0] - widget.canvasx(0)
                toprel = bbox[1] - widget.canvasy(0)
                bottomrel = bbox[3] - widget.canvasy(0)
	    else:
		# The widget is a text widget.  Place balloon under
                # the character closest to the mouse.  The positions
                # returned by bbox are relative to the text widget
                # window (ie the visible part of the text only).
		leftrel = bbox[0]
                toprel = bbox[1]
		bottomrel = bbox[1] + bbox[3]
	else:
	    leftrel = 0
            toprel = 0
	    bottomrel = widget.winfo_height()

        xpointer, ypointer = widget.winfo_pointerxy()   # -1 if off screen

        if xpointer >= 0 and self['relmouse'] in ('both', 'x'):
            x = xpointer
        else:
            x = leftrel + widget.winfo_rootx()
        x = x + self['xoffset']

        if ypointer >= 0 and self['relmouse'] in ('both', 'y'):
            y = ypointer
        else:
            y = bottomrel + widget.winfo_rooty()
        y = y + self['yoffset']

        edges = (string.atoi(str(self.cget('hull_highlightthickness'))) +
            string.atoi(str(self.cget('hull_borderwidth')))) * 2
        if x + self._label.winfo_reqwidth() + edges > screenWidth:
            x = screenWidth - self._label.winfo_reqwidth() - edges

        if y + self._label.winfo_reqheight() + edges > screenHeight:
            if ypointer >= 0 and self['relmouse'] in ('both', 'y'):
                y = ypointer
            else:
                y = toprel + widget.winfo_rooty()
            y = y - self._label.winfo_reqheight() - self['yoffset'] - edges

        setgeometryanddeiconify(self, '+%d+%d' % (x, y))

######################################################################
### File: PmwButtonBox.py
# Based on iwidgets2.2.0/buttonbox.itk code.

import types
import Tkinter


class ButtonBox(MegaWidget):
    def __init__(self, parent = None, **kw):

	# Define the megawidget options.
	
	optiondefs = (
	    ('labelmargin',       0,              INITOPT),
	    ('labelpos',          None,           INITOPT),
	    ('orient',            'horizontal',   INITOPT),
	    ('padx',              3,              INITOPT),
	    ('pady',              3,              INITOPT),
	)
	self.defineoptions(kw, optiondefs, dynamicGroups = ('Button',))

	# Initialise the base class (after defining the options).
	MegaWidget.__init__(self, parent)

	# Create the components.
	interior = self.interior()
	if self['labelpos'] is None:
	    self._buttonBoxFrame = self._hull
	    columnOrRow = 0
	else:
	    self._buttonBoxFrame = self.createcomponent('frame',
		    (), None,
		    Tkinter.Frame, (interior,))
	    self._buttonBoxFrame.grid(column=2, row=2, sticky='nsew')
	    columnOrRow = 2

	    self.createlabel(interior)

	orient = self['orient']
	if orient == 'horizontal':
	    interior.grid_columnconfigure(columnOrRow, weight = 1)
	elif orient == 'vertical':
	    interior.grid_rowconfigure(columnOrRow, weight = 1)
	else:
	    raise ValueError("bad orient option %s: must be either 'horizontal' or "
                             "'vertical'" % repr(orient))

	# Initialise instance variables.

	# List of tuples describing the buttons:
	#   - name
	#   - button widget
	self._buttonList = []

	# The index of the default button.
	self._defaultButton = None

	self._timerId = None

	# Check keywords and initialise options.
	self.initialiseoptions()

    def destroy(self):
	if self._timerId:
	    self.after_cancel(self._timerId)
	    self._timerId = None
	MegaWidget.destroy(self)

    def numbuttons(self):
        return len(self._buttonList)

    def index(self, index, forInsert = 0):
	listLength = len(self._buttonList)
        if isinstance(index, int):
	    if forInsert and index <= listLength:
		return index
	    elif not forInsert and index < listLength:
		return index
	    else:
		raise ValueError('index "%s" is out of range' % index)
	elif index is END:
	    if forInsert:
		return listLength
	    elif listLength > 0:
		return listLength - 1
	    else:
		raise ValueError('ButtonBox has no buttons')
	elif index is DEFAULT:
	    if self._defaultButton is not None:
		return self._defaultButton
	    raise ValueError('ButtonBox has no default')
	else:
            names = map(lambda t: t[0], self._buttonList)
            if index in names:
                return names.index(index)
	    validValues = 'a name, a number, END or DEFAULT'
	    raise ValueError('bad index "%s": must be %s' % (index, validValues))

    def insert(self, componentName, beforeComponent = 0, **kw):
	if componentName in self.components():
	    raise ValueError('button "%s" already exists' % componentName)
	if not kw.has_key('text'):
	    kw['text'] = componentName
        kw['default'] = 'normal'
	button = apply(self.createcomponent, (componentName,
		(), 'Button',
		Tkinter.Button, (self._buttonBoxFrame,)), kw)

	index = self.index(beforeComponent, 1)
	horizontal = self['orient'] == 'horizontal'
	numButtons = len(self._buttonList)

	# Shift buttons up one position.
	for i in range(numButtons - 1, index - 1, -1):
	    widget = self._buttonList[i][1]
	    pos = i * 2 + 3
	    if horizontal:
		widget.grid(column = pos, row = 0)
	    else:
		widget.grid(column = 0, row = pos)

	# Display the new button.
	if horizontal:
	    button.grid(column = index * 2 + 1, row = 0, sticky = 'ew',
		    padx = self['padx'], pady = self['pady'])
	    self._buttonBoxFrame.grid_columnconfigure(
		    numButtons * 2 + 2, weight = 1)
	else:
	    button.grid(column = 0, row = index * 2 + 1, sticky = 'ew',
		    padx = self['padx'], pady = self['pady'])
	    self._buttonBoxFrame.grid_rowconfigure(
		    numButtons * 2 + 2, weight = 1)
	self._buttonList.insert(index, (componentName, button))

	return button

    def add(self, componentName, **kw):
        return apply(self.insert, (componentName, len(self._buttonList)), kw)

    def delete(self, index):
        index = self.index(index)
	(name, widget) = self._buttonList[index]
	widget.grid_forget()
	self.destroycomponent(name)

	numButtons = len(self._buttonList)

	# Shift buttons down one position.
	horizontal = self['orient'] == 'horizontal'
	for i in range(index + 1, numButtons):
	    widget = self._buttonList[i][1]
	    pos = i * 2 - 1
	    if horizontal:
		widget.grid(column = pos, row = 0)
	    else:
		widget.grid(column = 0, row = pos)

	if horizontal:
	    self._buttonBoxFrame.grid_columnconfigure(numButtons * 2 - 1,
		    minsize = 0)
	    self._buttonBoxFrame.grid_columnconfigure(numButtons * 2, weight = 0)
	else:
	    self._buttonBoxFrame.grid_rowconfigure(numButtons * 2, weight = 0)
	del self._buttonList[index]

    def setdefault(self, index):
	# Turn off the default ring around the current default button.
	if self._defaultButton is not None:
	    button = self._buttonList[self._defaultButton][1]
	    button.configure(default = 'normal')
	    self._defaultButton = None

	# Turn on the default ring around the new default button.
	if index is not None:
	    index = self.index(index)
	    self._defaultButton = index
	    button = self._buttonList[index][1]
	    button.configure(default = 'active')

    def invoke(self, index = DEFAULT, noFlash = 0):
	# Invoke the callback associated with the *index* button.  If
	# *noFlash* is not set, flash the button to indicate to the
	# user that something happened.

	button = self._buttonList[self.index(index)][1]
	if not noFlash:
	    state = button.cget('state')
	    relief = button.cget('relief')
	    button.configure(state = 'active', relief = 'sunken')
	    self.update_idletasks()
	    self.after(100)
	    button.configure(state = state, relief = relief)
	return button.invoke()

    def button(self, buttonIndex):
	return self._buttonList[self.index(buttonIndex)][1]

    def alignbuttons(self, when = 'later'):
	if when == 'later':
	    if not self._timerId:
		self._timerId = self.after_idle(self.alignbuttons, 'now')
	    return
	self.update_idletasks()
	self._timerId = None

        # Determine the width of the maximum length button.
        max = 0
        horizontal = (self['orient'] == 'horizontal')
        for index in range(len(self._buttonList)):
            gridIndex = index * 2 + 1
            if horizontal:
                width = self._buttonBoxFrame.grid_bbox(gridIndex, 0)[2]
            else:
                width = self._buttonBoxFrame.grid_bbox(0, gridIndex)[2]
            if width > max:
                max = width

        # Set the width of all the buttons to be the same.
	if horizontal:
	    for index in range(len(self._buttonList)):
		self._buttonBoxFrame.grid_columnconfigure(index * 2 + 1,
			minsize = max)
	else:
	    self._buttonBoxFrame.grid_columnconfigure(0, minsize = max)

######################################################################
### File: PmwEntryField.py
# Based on iwidgets2.2.0/entryfield.itk code.

import re
import string
import types
import Tkinter


# Possible return values of validation functions.
OK = 1
ERROR = 0
PARTIAL = -1

class EntryField(MegaWidget):
    _classBindingsDefinedFor = 0

    def __init__(self, parent = None, **kw):

	# Define the megawidget options.
	
	optiondefs = (
	    ('command',           None,        None),
	    ('errorbackground',   'pink',      None),
	    ('invalidcommand',    self.bell,   None),
	    ('labelmargin',       0,           INITOPT),
	    ('labelpos',          None,        INITOPT),
	    ('modifiedcommand',   None,        None),
	    ('sticky',            'ew',        INITOPT),
	    ('validate',          None,        self._validate),
	    ('extravalidators',   {},          None),
	    ('value',             '',          INITOPT),
	)
	self.defineoptions(kw, optiondefs)

	# Initialise the base class (after defining the options).
	MegaWidget.__init__(self, parent)

	# Create the components.
	interior = self.interior()
	self._entryFieldEntry = self.createcomponent('entry',
		(), None,
		Tkinter.Entry, (interior,))
	self._entryFieldEntry.grid(column=2, row=2, sticky=self['sticky'])
	if self['value'] != '':
	    self.__setEntry(self['value'])
	interior.grid_columnconfigure(2, weight=1)
	interior.grid_rowconfigure(2, weight=1)

	self.createlabel(interior)

	# Initialise instance variables.

	self.normalBackground = None
        self._previousText = None

	# Initialise instance.

	_registerEntryField(self._entryFieldEntry, self)

        # Establish the special class bindings if not already done.
        # Also create bindings if the Tkinter default interpreter has
        # changed.  Use Tkinter._default_root to create class
        # bindings, so that a reference to root is created by
        # bind_class rather than a reference to self, which would
        # prevent object cleanup.
        if EntryField._classBindingsDefinedFor != Tkinter._default_root:
	    tagList = self._entryFieldEntry.bindtags()
            root  = Tkinter._default_root
	    	    
	    allSequences = {}
	    for tag in tagList:

                sequences = root.bind_class(tag)
                if isinstance(sequences, basestring):
                    # In old versions of Tkinter, bind_class returns a string
                    sequences = root.tk.splitlist(sequences)

		for sequence in sequences:
		    allSequences[sequence] = None
	    for sequence in allSequences.keys():
		root.bind_class('EntryFieldPre', sequence, _preProcess)
		root.bind_class('EntryFieldPost', sequence, _postProcess)

	    EntryField._classBindingsDefinedFor = root

	self._entryFieldEntry.bindtags(('EntryFieldPre',) +
		self._entryFieldEntry.bindtags() + ('EntryFieldPost',))
	self._entryFieldEntry.bind('<Return>', self._executeCommand)

	# Check keywords and initialise options.
	self.initialiseoptions()

    def destroy(self):
	_deregisterEntryField(self._entryFieldEntry)
        MegaWidget.destroy(self)

    def _getValidatorFunc(self, validator, index):
	# Search the extra and standard validator lists for the
	# given 'validator'.  If 'validator' is an alias, then
	# continue the search using the alias.  Make sure that
	# self-referencial aliases do not cause infinite loops.

	extraValidators = self['extravalidators']
	traversedValidators = []

	while 1:
	    traversedValidators.append(validator)
	    if extraValidators.has_key(validator):
		validator = extraValidators[validator][index]
	    elif _standardValidators.has_key(validator):
		validator = _standardValidators[validator][index]
	    else:
		return validator
	    if validator in traversedValidators:
		return validator

    def _validate(self):
	dict = {
	    'validator' : None,
	    'min' : None,
	    'max' : None,
	    'minstrict' : 1,
	    'maxstrict' : 1,
	}
	opt = self['validate']
        if isinstance(opt, dict):
	    dict.update(opt)
	else:
	    dict['validator'] = opt

	# Look up validator maps and replace 'validator' field with
	# the corresponding function.
	validator = dict['validator']
	valFunction = self._getValidatorFunc(validator, 0)
	self._checkValidateFunction(valFunction, 'validate', validator)
	dict['validator'] = valFunction

	# Look up validator maps and replace 'stringtovalue' field
	# with the corresponding function.
	if dict.has_key('stringtovalue'):
	    stringtovalue = dict['stringtovalue'] 
	    strFunction = self._getValidatorFunc(stringtovalue, 1)
	    self._checkValidateFunction(
		    strFunction, 'stringtovalue', stringtovalue)
	else:
	    strFunction = self._getValidatorFunc(validator, 1)
	    if strFunction == validator:
		strFunction = len
	dict['stringtovalue'] = strFunction

	self._validationInfo = dict
	args = dict.copy()
	del args['validator']
	del args['min']
	del args['max']
	del args['minstrict']
	del args['maxstrict']
	del args['stringtovalue']
	self._validationArgs = args
        self._previousText = None

        if isinstance(dict['min'], basestring) and strFunction is not None:
	    dict['min'] = apply(strFunction, (dict['min'],), args)
        if isinstance(dict['max'], basestring) and strFunction is not None:
	    dict['max'] = apply(strFunction, (dict['max'],), args)

	self._checkValidity()

    def _checkValidateFunction(self, function, option, validator):
	# Raise an error if 'function' is not a function or None.

	if function is not None and not callable(function):
	    extraValidators = self['extravalidators']
	    extra = extraValidators.keys()
	    extra.sort()
	    extra = tuple(extra)
	    standard = _standardValidators.keys()
	    standard.sort()
	    standard = tuple(standard)
	    msg = 'bad %s value "%s":  must be a function or one of ' \
		'the standard validators %s or extra validators %s'
	    raise ValueError(msg % (option, validator, standard, extra))

    def _executeCommand(self, event = None):
	cmd = self['command']
	if callable(cmd):
            if event is None:
                # Return result of command for invoke() method.
                return cmd()
            else:
                cmd()
	    
    def _preProcess(self):

        self._previousText = self._entryFieldEntry.get()
        self._previousICursor = self._entryFieldEntry.index('insert')
        self._previousXview = self._entryFieldEntry.index('@0')
	if self._entryFieldEntry.selection_present():
	    self._previousSel= (self._entryFieldEntry.index('sel.first'),
		self._entryFieldEntry.index('sel.last'))
	else:
	    self._previousSel = None

    def _postProcess(self):

	# No need to check if text has not changed.
	previousText = self._previousText
	if previousText == self._entryFieldEntry.get():
	    return self.valid()

	valid = self._checkValidity()
        if self.hulldestroyed():
            # The invalidcommand called by _checkValidity() destroyed us.
            return valid

	cmd = self['modifiedcommand']
	if callable(cmd) and previousText != self._entryFieldEntry.get():
	    cmd()
	return valid
	    
    def checkentry(self):
	# If there is a variable specified by the entry_textvariable
	# option, checkentry() should be called after the set() method
	# of the variable is called.

	self._previousText = None
	return self._postProcess()

    def _getValidity(self):
	text = self._entryFieldEntry.get()
	dict = self._validationInfo
	args = self._validationArgs

	if dict['validator'] is not None:
	    status = apply(dict['validator'], (text,), args)
	    if status != OK:
		return status

	# Check for out of (min, max) range.
	if dict['stringtovalue'] is not None:
	    min = dict['min']
	    max = dict['max']
	    if min is None and max is None:
		return OK
	    val = apply(dict['stringtovalue'], (text,), args)
	    if min is not None and val < min:
		if dict['minstrict']:
		    return ERROR
		else:
		    return PARTIAL
	    if max is not None and val > max:
		if dict['maxstrict']:
		    return ERROR
		else:
		    return PARTIAL
	return OK

    def _checkValidity(self):
	valid = self._getValidity()
	oldValidity = valid

	if valid == ERROR:
	    # The entry is invalid.
	    cmd = self['invalidcommand']
	    if callable(cmd):
		cmd()
            if self.hulldestroyed():
                # The invalidcommand destroyed us.
                return oldValidity

	    # Restore the entry to its previous value.
	    if self._previousText is not None:
		self.__setEntry(self._previousText)
		self._entryFieldEntry.icursor(self._previousICursor)
		self._entryFieldEntry.xview(self._previousXview)
		if self._previousSel is not None:
		    self._entryFieldEntry.selection_range(self._previousSel[0],
			self._previousSel[1])

		# Check if the saved text is valid as well.
		valid = self._getValidity()

	self._valid = valid

        if self.hulldestroyed():
            # The validator or stringtovalue commands called by
            # _checkValidity() destroyed us.
            return oldValidity

	if valid == OK:
	    if self.normalBackground is not None:
		self._entryFieldEntry.configure(
			background = self.normalBackground)
		self.normalBackground = None
	else:
	    if self.normalBackground is None:
		self.normalBackground = self._entryFieldEntry.cget('background')
		self._entryFieldEntry.configure(
			background = self['errorbackground'])

        return oldValidity

    def invoke(self):
	return self._executeCommand()

    def valid(self):
        return self._valid == OK

    def clear(self):
        self.setentry('')

    def __setEntry(self, text):
	oldState = str(self._entryFieldEntry.cget('state'))
	if oldState != 'normal':
	    self._entryFieldEntry.configure(state='normal')
	self._entryFieldEntry.delete(0, 'end')
	self._entryFieldEntry.insert(0, text)
	if oldState != 'normal':
	    self._entryFieldEntry.configure(state=oldState)

    def setentry(self, text):
	self._preProcess()
        self.__setEntry(text)
	return self._postProcess()

    def getvalue(self):
        return self._entryFieldEntry.get()

    def setvalue(self, text):
        return self.setentry(text)

forwardmethods(EntryField, Tkinter.Entry, '_entryFieldEntry')

# ======================================================================


# Entry field validation functions

_numericregex = re.compile('^[0-9]*$')
_alphabeticregex = re.compile('^[a-z]*$', re.IGNORECASE)
_alphanumericregex = re.compile('^[0-9a-z]*$', re.IGNORECASE)

def numericvalidator(text):
    if text == '':
        return PARTIAL
    else:
	if _numericregex.match(text) is None:
	    return ERROR
	else:
	    return OK
    
def integervalidator(text):
    if text in ('', '-', '+'):
        return PARTIAL
    try:
	string.atol(text)
	return OK
    except ValueError:
	return ERROR
    
def alphabeticvalidator(text):
    if _alphabeticregex.match(text) is None:
	return ERROR
    else:
	return OK
    
def alphanumericvalidator(text):
    if _alphanumericregex.match(text) is None:
	return ERROR
    else:
	return OK
    
def hexadecimalvalidator(text):
    if text in ('', '0x', '0X', '+', '+0x', '+0X', '-', '-0x', '-0X'):
        return PARTIAL
    try:
	string.atol(text, 16)
	return OK
    except ValueError:
	return ERROR
    
def realvalidator(text, separator = '.'):
    if separator != '.':
	if string.find(text, '.') >= 0:
	    return ERROR
	index = string.find(text, separator)
	if index >= 0:
	    text = text[:index] + '.' + text[index + 1:]
    try:
	string.atof(text)
	return OK
    except ValueError:
	# Check if the string could be made valid by appending a digit
	# eg ('-', '+', '.', '-.', '+.', '1.23e', '1E-').
	if len(text) == 0:
	    return PARTIAL
	if text[-1] in string.digits:
	    return ERROR
	try:
	    string.atof(text + '0')
	    return PARTIAL
	except ValueError:
	    return ERROR
    
def timevalidator(text, separator = ':'):
    try:
	timestringtoseconds(text, separator)
	return OK
    except ValueError:
	if len(text) > 0 and text[0] in ('+', '-'):
	    text = text[1:]
	if re.search('[^0-9' + separator + ']', text) is not None:
	    return ERROR
	return PARTIAL

def datevalidator(text, format = 'ymd', separator = '/'):
    try:
	datestringtojdn(text, format, separator)
	return OK
    except ValueError:
	if re.search('[^0-9' + separator + ']', text) is not None:
	    return ERROR
	return PARTIAL

_standardValidators = {
    'numeric'      : (numericvalidator,      string.atol),
    'integer'      : (integervalidator,      string.atol),
    'hexadecimal'  : (hexadecimalvalidator,  lambda s: string.atol(s, 16)),
    'real'         : (realvalidator,         stringtoreal),
    'alphabetic'   : (alphabeticvalidator,   len),
    'alphanumeric' : (alphanumericvalidator, len),
    'time'         : (timevalidator,         timestringtoseconds),
    'date'         : (datevalidator,         datestringtojdn),
}

_entryCache = {}

def _registerEntryField(entry, entryField):
    # Register an EntryField widget for an Entry widget

    _entryCache[entry] = entryField

def _deregisterEntryField(entry):
    # Deregister an Entry widget
    del _entryCache[entry]

def _preProcess(event):
    # Forward preprocess events for an Entry to it's EntryField

    _entryCache[event.widget]._preProcess()

def _postProcess(event):
    # Forward postprocess events for an Entry to it's EntryField

    # The function specified by the 'command' option may have destroyed
    # the megawidget in a binding earlier in bindtags, so need to check.
    if _entryCache.has_key(event.widget):
        _entryCache[event.widget]._postProcess()

######################################################################
### File: PmwGroup.py
import string
import Tkinter


def aligngrouptags(groups):
    # Adjust the y position of the tags in /groups/ so that they all
    # have the height of the highest tag.

    maxTagHeight = 0
    for group in groups:
	if group._tag is None:
	    height = (string.atoi(str(group._ring.cget('borderwidth'))) +
		    string.atoi(str(group._ring.cget('highlightthickness'))))
	else:
	    height = group._tag.winfo_reqheight()
	if maxTagHeight < height:
	    maxTagHeight = height

    for group in groups:
	ringBorder = (string.atoi(str(group._ring.cget('borderwidth'))) +
		string.atoi(str(group._ring.cget('highlightthickness'))))
	topBorder = maxTagHeight / 2 - ringBorder / 2
	group._hull.grid_rowconfigure(0, minsize = topBorder)
	group._ring.grid_rowconfigure(0,
		minsize = maxTagHeight - topBorder - ringBorder)
	if group._tag is not None:
	    group._tag.place(y = maxTagHeight / 2)

class Group( MegaWidget ):
    def __init__(self, parent = None, **kw):

        # Define the megawidget options.
	
        optiondefs = (
	    ('collapsedsize',    6,         INITOPT),
	    ('ring_borderwidth', 2,         None),
	    ('ring_relief',      'groove',  None),
	    ('tagindent',        10,        INITOPT),
        )
        self.defineoptions(kw, optiondefs)

        # Initialise the base class (after defining the options).
        MegaWidget.__init__(self, parent)

        # Create the components.
        interior = MegaWidget.interior(self)

	self._ring = self.createcomponent(
	    'ring', 
	    (), None,
	    Tkinter.Frame, (interior,), 
	    )

	self._groupChildSite = self.createcomponent(
	    'groupchildsite',
	    (), None,
	    Tkinter.Frame, (self._ring,)
	    )

        self._tag = self.createcomponent(
	    'tag',
	    (), None,
	    Tkinter.Label, (interior,),
	    )

	ringBorder = (string.atoi(str(self._ring.cget('borderwidth'))) +
		string.atoi(str(self._ring.cget('highlightthickness'))))
	if self._tag is None:
	    tagHeight = ringBorder
	else:
	    tagHeight = self._tag.winfo_reqheight()
	    self._tag.place(
		    x = ringBorder + self['tagindent'],
		    y = tagHeight / 2,
		    anchor = 'w')

	topBorder = tagHeight / 2 - ringBorder / 2
	self._ring.grid(column = 0, row = 1, sticky = 'nsew')
	interior.grid_columnconfigure(0, weight = 1)
	interior.grid_rowconfigure(1, weight = 1)
	interior.grid_rowconfigure(0, minsize = topBorder)

	self._groupChildSite.grid(column = 0, row = 1, sticky = 'nsew')
	self._ring.grid_columnconfigure(0, weight = 1)
	self._ring.grid_rowconfigure(1, weight = 1)
	self._ring.grid_rowconfigure(0,
		minsize = tagHeight - topBorder - ringBorder)

        self.showing = 1

        # Check keywords and initialise options.
        self.initialiseoptions()

    def toggle(self):
        if self.showing:
            self.collapse()
        else:
            self.expand()
        self.showing = not self.showing

    def expand(self):
        self._groupChildSite.grid(column = 0, row = 1, sticky = 'nsew')

    def collapse(self):
        self._groupChildSite.grid_forget()
	if self._tag is None:
	    tagHeight = 0
	else:
            tagHeight = self._tag.winfo_reqheight()
        self._ring.configure(height=(tagHeight / 2) + self['collapsedsize'])

    def interior(self):
        return self._groupChildSite

######################################################################
### File: PmwLabeledWidget.py
import Tkinter


class LabeledWidget(MegaWidget):
    def __init__(self, parent = None, **kw):

	# Define the megawidget options.
	
	optiondefs = (
	    ('labelmargin',            0,      INITOPT),
	    ('labelpos',               None,   INITOPT),
	    ('sticky',                 'nsew', INITOPT),
	)
	self.defineoptions(kw, optiondefs)

	# Initialise the base class (after defining the options).
	MegaWidget.__init__(self, parent)

	# Create the components.
	interior = MegaWidget.interior(self)
	self._labelChildSite = self.createcomponent('labelchildsite',
		(), None,
		Tkinter.Frame, (interior,))
	self._labelChildSite.grid(column=2, row=2, sticky=self['sticky'])
	interior.grid_columnconfigure(2, weight=1)
	interior.grid_rowconfigure(2, weight=1)

	self.createlabel(interior)

	# Check keywords and initialise options.
	self.initialiseoptions()

    def interior(self):
	return self._labelChildSite

######################################################################
### File: PmwMainMenuBar.py
# Main menubar

import string
import types
import Tkinter


class MainMenuBar(MegaArchetype):

    def __init__(self, parent = None, **kw):

        # Define the megawidget options.
        
        optiondefs = (
            ('balloon',      None,       None),
            ('hotkeys',      1,          INITOPT),
            ('hull_tearoff', 0,          None),
        )
        self.defineoptions(kw, optiondefs, dynamicGroups = ('Menu',))

        # Initialise the base class (after defining the options).
        MegaArchetype.__init__(self, parent, Tkinter.Menu)

        self._menuInfo = {}
        self._menuInfo[None] = (None, [])
        # Map from a menu name to a tuple of information about the menu.
        # The first item in the tuple is the name of the parent menu (for
        # toplevel menus this is None). The second item in the tuple is
        # a list of status help messages for each item in the menu.
        # The key for the information for the main menubar is None.

        self._menu = self.interior()
        self._menu.bind('<Leave>', self._resetHelpmessage)
        self._menu.bind('<Motion>',
            lambda event=None, self=self: self._menuHelp(event, None))

        # Check keywords and initialise options.
        self.initialiseoptions()

    def deletemenuitems(self, menuName, start, end = None):
        self.component(menuName).delete(start, end)
        if end is None:
            del self._menuInfo[menuName][1][start]
        else:
            self._menuInfo[menuName][1][start:end+1] = []

    def deletemenu(self, menuName):
        """Delete should be called for cascaded menus before main menus.
        """

        parentName = self._menuInfo[menuName][0]
        del self._menuInfo[menuName]
        if parentName is None:
            parentMenu = self._menu
        else:
            parentMenu = self.component(parentName)

        menu = self.component(menuName)
        menuId = str(menu)
        for item in range(parentMenu.index('end') + 1):
            if parentMenu.type(item) == 'cascade':
                itemMenu = str(parentMenu.entrycget(item, 'menu'))
                if itemMenu == menuId:
                    parentMenu.delete(item)
                    del self._menuInfo[parentName][1][item]
                    break

        self.destroycomponent(menuName)

    def disableall(self):
        for index in range(len(self._menuInfo[None][1])):
            self.entryconfigure(index, state = 'disabled')

    def enableall(self):
        for index in range(len(self._menuInfo[None][1])):
            self.entryconfigure(index, state = 'normal')

    def addmenu(self, menuName, balloonHelp, statusHelp = None,
            traverseSpec = None, **kw):
        if statusHelp is None:
            statusHelp = balloonHelp
        self._addmenu(None, menuName, balloonHelp, statusHelp,
            traverseSpec, kw)

    def addcascademenu(self, parentMenuName, menuName, statusHelp='',
            traverseSpec = None, **kw):
        self._addmenu(parentMenuName, menuName, None, statusHelp,
            traverseSpec, kw)

    def _addmenu(self, parentMenuName, menuName, balloonHelp, statusHelp,
            traverseSpec, kw):

        if (menuName) in self.components():
            raise ValueError('menu "%s" already exists' % menuName)

        menukw = {}
        if kw.has_key('tearoff'):
            menukw['tearoff'] = kw['tearoff']
            del kw['tearoff']
        else:
            menukw['tearoff'] = 0
        if kw.has_key('name'):
            menukw['name'] = kw['name']
            del kw['name']

        if not kw.has_key('label'):
            kw['label'] = menuName

        self._addHotkeyToOptions(parentMenuName, kw, traverseSpec)

        if parentMenuName is None:
            parentMenu = self._menu
            balloon = self['balloon']
            # Bug in Tk: balloon help not implemented
            # if balloon is not None:
            #     balloon.mainmenubind(parentMenu, balloonHelp, statusHelp)
        else:
            parentMenu = self.component(parentMenuName)

        apply(parentMenu.add_cascade, (), kw)

        menu = apply(self.createcomponent, (menuName,
                (), 'Menu',
                Tkinter.Menu, (parentMenu,)), menukw)
        parentMenu.entryconfigure('end', menu = menu)

        self._menuInfo[parentMenuName][1].append(statusHelp)
        self._menuInfo[menuName] = (parentMenuName, [])

        menu.bind('<Leave>', self._resetHelpmessage)
        menu.bind('<Motion>', 
            lambda event=None, self=self, menuName=menuName:
                    self._menuHelp(event, menuName))

    def addmenuitem(self, menuName, itemType, statusHelp = '',
            traverseSpec = None, **kw):

        menu = self.component(menuName)
        if itemType != 'separator':
            self._addHotkeyToOptions(menuName, kw, traverseSpec)

        if itemType == 'command':
            command = menu.add_command
        elif itemType == 'separator':
            command = menu.add_separator
        elif itemType == 'checkbutton':
            command = menu.add_checkbutton
        elif itemType == 'radiobutton':
            command = menu.add_radiobutton
        elif itemType == 'cascade':
            command = menu.add_cascade
        else:
            raise ValueError('unknown menuitem type "%s"' % itemType)

        self._menuInfo[menuName][1].append(statusHelp)
        apply(command, (), kw)

    def _addHotkeyToOptions(self, menuName, kw, traverseSpec):

        if (not self['hotkeys'] or kw.has_key('underline') or
                not kw.has_key('label')):
            return

        if isinstance(traverseSpec, int):
            kw['underline'] = traverseSpec
            return

        if menuName is None:
            menu = self._menu
        else:
            menu = self.component(menuName)
        hotkeyList = []
        end = menu.index('end')
        if end is not None:
            for item in range(end + 1):
                if menu.type(item) not in ('separator', 'tearoff'):
                    underline = \
                            string.atoi(str(menu.entrycget(item, 'underline')))
                    if underline != -1:
                        label = str(menu.entrycget(item, 'label'))
                        if underline < len(label):
                            hotkey = string.lower(label[underline])
                            if hotkey not in hotkeyList:
                                hotkeyList.append(hotkey)

        name = kw['label']

        if isinstance(traverseSpec, basestring):
            lowerLetter = string.lower(traverseSpec)
            if traverseSpec in name and lowerLetter not in hotkeyList:
                kw['underline'] = string.index(name, traverseSpec)
        else:
            targets = string.digits + string.letters
            lowerName = string.lower(name)
            for letter_index in range(len(name)):
                letter = lowerName[letter_index]
                if letter in targets and letter not in hotkeyList:
                    kw['underline'] = letter_index
                    break

    def _menuHelp(self, event, menuName):
        if menuName is None:
            menu = self._menu
            index = menu.index('@%d'% event.x)
        else:
            menu = self.component(menuName)
            index = menu.index('@%d'% event.y)

        balloon = self['balloon']
        if balloon is not None:
            if index is None:
                balloon.showstatus('')
            else:
                if str(menu.cget('tearoff')) == '1':
                    index = index - 1
                if index >= 0:
                    help = self._menuInfo[menuName][1][index]
                    balloon.showstatus(help)

    def _resetHelpmessage(self, event=None):
        balloon = self['balloon']
        if balloon is not None:
            balloon.clearstatus()

forwardmethods(MainMenuBar, Tkinter.Menu, '_hull')

######################################################################
### File: PmwMenuBar.py
# Manager widget for menus.

import string
import types
import Tkinter


class MenuBar(MegaWidget):

    def __init__(self, parent = None, **kw):

        # Define the megawidget options.
        
        optiondefs = (
            ('balloon',      None,       None),
            ('hotkeys',      1,          INITOPT),
            ('padx',         0,          INITOPT),
        )
        self.defineoptions(kw, optiondefs, dynamicGroups = ('Menu', 'Button'))

        # Initialise the base class (after defining the options).
        MegaWidget.__init__(self, parent)

        self._menuInfo = {}
        # Map from a menu name to a tuple of information about the menu.
        # The first item in the tuple is the name of the parent menu (for
        # toplevel menus this is None). The second item in the tuple is
        # a list of status help messages for each item in the menu.
        # The third item in the tuple is the id of the binding used
        # to detect mouse motion to display status help.
        # Information for the toplevel menubuttons is not stored here.

        self._mydeletecommand = self.component('hull').tk.deletecommand
        # Cache this method for use later.

        # Check keywords and initialise options.
        self.initialiseoptions()

    def deletemenuitems(self, menuName, start, end = None):
        self.component(menuName + '-menu').delete(start, end)
        if end is None:
            del self._menuInfo[menuName][1][start]
        else:
            self._menuInfo[menuName][1][start:end+1] = []

    def deletemenu(self, menuName):
        """Delete should be called for cascaded menus before main menus.
        """

        # Clean up binding for this menu.
        parentName = self._menuInfo[menuName][0]
        bindId = self._menuInfo[menuName][2]
        _bindtag = 'PmwMenuBar' + str(self) + menuName
        self.unbind_class(_bindtag, '<Motion>')
        self._mydeletecommand(bindId) # unbind_class does not clean up
        del self._menuInfo[menuName]

        if parentName is None:
            self.destroycomponent(menuName + '-button')
        else:
            parentMenu = self.component(parentName + '-menu')

            menu = self.component(menuName + '-menu')
            menuId = str(menu)
            for item in range(parentMenu.index('end') + 1):
                if parentMenu.type(item) == 'cascade':
                    itemMenu = str(parentMenu.entrycget(item, 'menu'))
                    if itemMenu == menuId:
                        parentMenu.delete(item)
                        del self._menuInfo[parentName][1][item]
                        break

        self.destroycomponent(menuName + '-menu')

    def disableall(self):
        for menuName in self._menuInfo.keys():
            if self._menuInfo[menuName][0] is None:
                menubutton = self.component(menuName + '-button')
                menubutton.configure(state = 'disabled')

    def enableall(self):
        for menuName in self._menuInfo.keys():
            if self._menuInfo[menuName][0] is None:
                menubutton = self.component(menuName + '-button')
                menubutton.configure(state = 'normal')

    def addmenu(self, menuName, balloonHelp, statusHelp = None,
            side = 'left', traverseSpec = None, **kw):

        self._addmenu(None, menuName, balloonHelp, statusHelp,
            traverseSpec, side, 'text', kw)

    def addcascademenu(self, parentMenuName, menuName, statusHelp = '',
            traverseSpec = None, **kw):

        self._addmenu(parentMenuName, menuName, None, statusHelp,
            traverseSpec, None, 'label', kw)

    def _addmenu(self, parentMenuName, menuName, balloonHelp, statusHelp,
            traverseSpec, side, textKey, kw):

        if (menuName + '-menu') in self.components():
            raise ValueError('menu "%s" already exists' % menuName)

        menukw = {}
        if kw.has_key('tearoff'):
            menukw['tearoff'] = kw['tearoff']
            del kw['tearoff']
        else:
            menukw['tearoff'] = 0

        if not kw.has_key(textKey):
            kw[textKey] = menuName

        self._addHotkeyToOptions(parentMenuName, kw, textKey, traverseSpec)

        if parentMenuName is None:
            button = apply(self.createcomponent, (menuName + '-button',
                    (), 'Button',
                    Tkinter.Menubutton, (self.interior(),)), kw)
            button.pack(side=side, padx = self['padx'])
            balloon = self['balloon']
            if balloon is not None:
                balloon.bind(button, balloonHelp, statusHelp)
            parentMenu = button
        else:
            parentMenu = self.component(parentMenuName + '-menu')
            apply(parentMenu.add_cascade, (), kw)
            self._menuInfo[parentMenuName][1].append(statusHelp)

        menu = apply(self.createcomponent, (menuName + '-menu',
                (), 'Menu',
                Tkinter.Menu, (parentMenu,)), menukw)
        if parentMenuName is None:
            button.configure(menu = menu)
        else:
            parentMenu.entryconfigure('end', menu = menu)

        # Need to put this binding after the class bindings so that
        # menu.index() does not lag behind.
        _bindtag = 'PmwMenuBar' + str(self) + menuName
        bindId = self.bind_class(_bindtag, '<Motion>',
            lambda event=None, self=self, menuName=menuName:
                    self._menuHelp(menuName))
        menu.bindtags(menu.bindtags() + (_bindtag,))
        menu.bind('<Leave>', self._resetHelpmessage)

        self._menuInfo[menuName] = (parentMenuName, [], bindId)

    def addmenuitem(self, menuName, itemType, statusHelp = '',
            traverseSpec = None, **kw):

        menu = self.component(menuName + '-menu')
        if itemType != 'separator':
            self._addHotkeyToOptions(menuName, kw, 'label', traverseSpec)

        if itemType == 'command':
            command = menu.add_command
        elif itemType == 'separator':
            command = menu.add_separator
        elif itemType == 'checkbutton':
            command = menu.add_checkbutton
        elif itemType == 'radiobutton':
            command = menu.add_radiobutton
        elif itemType == 'cascade':
            command = menu.add_cascade
        else:
            raise ValueError('unknown menuitem type "%s"' % itemType)

        self._menuInfo[menuName][1].append(statusHelp)
        apply(command, (), kw)

    def _addHotkeyToOptions(self, menuName, kw, textKey, traverseSpec):

        if (not self['hotkeys'] or kw.has_key('underline') or
                not kw.has_key(textKey)):
            return

        if isinstance(traverseSpec, int):
            kw['underline'] = traverseSpec
            return

        hotkeyList = []
        if menuName is None:
            for menuName in self._menuInfo.keys():
                if self._menuInfo[menuName][0] is None:
                    menubutton = self.component(menuName + '-button')
                    underline = string.atoi(str(menubutton.cget('underline')))
                    if underline != -1:
                        label = str(menubutton.cget(textKey))
                        if underline < len(label):
                            hotkey = string.lower(label[underline])
                            if hotkey not in hotkeyList:
                                hotkeyList.append(hotkey)
        else:
            menu = self.component(menuName + '-menu')
            end = menu.index('end')
            if end is not None:
                for item in range(end + 1):
                    if menu.type(item) not in ('separator', 'tearoff'):
                        underline = string.atoi(
                            str(menu.entrycget(item, 'underline')))
                        if underline != -1:
                            label = str(menu.entrycget(item, textKey))
                            if underline < len(label):
                                hotkey = string.lower(label[underline])
                                if hotkey not in hotkeyList:
                                    hotkeyList.append(hotkey)

        name = kw[textKey]

        if isinstance(traverseSpec, basestring):
            lowerLetter = string.lower(traverseSpec)
            if traverseSpec in name and lowerLetter not in hotkeyList:
                kw['underline'] = string.index(name, traverseSpec)
        else:
            targets = string.digits + string.letters
            lowerName = string.lower(name)
            for letter_index in range(len(name)):
                letter = lowerName[letter_index]
                if letter in targets and letter not in hotkeyList:
                    kw['underline'] = letter_index
                    break

    def _menuHelp(self, menuName):
        menu = self.component(menuName + '-menu')
        index = menu.index('active')

        balloon = self['balloon']
        if balloon is not None:
            if index is None:
                balloon.showstatus('')
            else:
                if str(menu.cget('tearoff')) == '1':
                    index = index - 1
                if index >= 0:
                    help = self._menuInfo[menuName][1][index]
                    balloon.showstatus(help)

    def _resetHelpmessage(self, event=None):
        balloon = self['balloon']
        if balloon is not None:
            balloon.clearstatus()

######################################################################
### File: PmwMessageBar.py
# Class to display messages in an information line.

import string
import Tkinter


class MessageBar(MegaWidget):
    def __init__(self, parent = None, **kw):

	# Define the megawidget options.
	
	defaultMessageTypes = {
			   # (priority, showtime, bells, logmessage)
	    'systemerror'  : (5, 10, 2, 1),
	    'usererror'    : (4, 5, 1, 0),
	    'busy'         : (3, 0, 0, 0),
	    'systemevent'  : (2, 5, 0, 0),
	    'userevent'    : (2, 5, 0, 0),
	    'help'         : (1, 5, 0, 0),
	    'state'        : (0, 0, 0, 0),
	}
	optiondefs = (
	    ('labelmargin',    0,                     INITOPT),
	    ('labelpos',       None,                  INITOPT),
	    ('messagetypes',   defaultMessageTypes,   INITOPT),
	    ('silent',         0,                     None),
	    ('sticky',         'ew',                  INITOPT),
	)
	self.defineoptions(kw, optiondefs)

	# Initialise the base class (after defining the options).
	MegaWidget.__init__(self, parent)

	# Create the components.
	interior = self.interior()
	self._messageBarEntry = self.createcomponent('entry',
		(), None,
		Tkinter.Entry, (interior,))

        # Can't always use 'disabled', since this greys out text in Tk 8.4.2
        try:
            self._messageBarEntry.configure(state = 'readonly')
        except Tkinter.TclError:
            self._messageBarEntry.configure(state = 'disabled')

	self._messageBarEntry.grid(column=2, row=2, sticky=self['sticky'])
	interior.grid_columnconfigure(2, weight=1)
	interior.grid_rowconfigure(2, weight=1)

	self.createlabel(interior)

	# Initialise instance variables.
	self._numPriorities = 0
	for info in self['messagetypes'].values():
	    if self._numPriorities < info[0]:
		self._numPriorities = info[0]

	self._numPriorities = self._numPriorities + 1
	self._timer = [None] * self._numPriorities
	self._messagetext = [''] * self._numPriorities
	self._activemessage = [0] * self._numPriorities

	# Check keywords and initialise options.
	self.initialiseoptions()

    def destroy(self):
	for timerId in self._timer:
	    if timerId is not None:
		self.after_cancel(timerId)
	self._timer = [None] * self._numPriorities
	MegaWidget.destroy(self)

    def message(self, type, text):
        # Display a message in the message bar.

	(priority, showtime, bells, logmessage) = self['messagetypes'][type]

	if not self['silent']:
	    for i in range(bells):
		if i != 0:
		    self.after(100)
		self.bell()

	self._activemessage[priority] = 1
	if text is None:
	    text = ''
	self._messagetext[priority] = string.replace(text, '\n', ' ')
	self._redisplayInfoMessage()

	if logmessage:
	    # Should log this text to a text widget.
	    pass

	if showtime > 0:
	    if self._timer[priority] is not None:
		self.after_cancel(self._timer[priority])

	    # Define a callback to clear this message after a time.
	    def _clearmessage(self=self, priority=priority):
		self._clearActivemessage(priority)

	    mseconds = int(showtime * 1000)
	    self._timer[priority] = self.after(mseconds, _clearmessage)

    def helpmessage(self, text):
        if text is None:
            self.resetmessages('help')
        else:
            self.message('help', text)

    def resetmessages(self, type):
	priority = self['messagetypes'][type][0]
	self._clearActivemessage(priority)
	for messagetype, info in self['messagetypes'].items():
	    thisPriority = info[0]
	    showtime = info[1]
	    if thisPriority < priority and showtime != 0:
		self._clearActivemessage(thisPriority)

    def _clearActivemessage(self, priority):
	self._activemessage[priority] = 0
	if self._timer[priority] is not None:
	    self.after_cancel(self._timer[priority])
	    self._timer[priority] = None
	self._redisplayInfoMessage()

    def _redisplayInfoMessage(self):
	text = ''
        for priority in range(self._numPriorities - 1, -1, -1):
	    if self._activemessage[priority]:
		text = self._messagetext[priority]
	        break
	self._messageBarEntry.configure(state = 'normal')
	self._messageBarEntry.delete(0, 'end')
	self._messageBarEntry.insert('end', text)

        # Can't always use 'disabled', since this greys out text in Tk 8.4.2
        try:
            self._messageBarEntry.configure(state = 'readonly')
        except Tkinter.TclError:
            self._messageBarEntry.configure(state = 'disabled')

forwardmethods(MessageBar, Tkinter.Entry, '_messageBarEntry')

######################################################################
### File: PmwMessageDialog.py
# Based on iwidgets2.2.0/messagedialog.itk code.

import Tkinter


class MessageDialog(Dialog):
    def __init__(self, parent = None, **kw):

	# Define the megawidget options.
	
	optiondefs = (
	    ('borderx',       20,    INITOPT),
	    ('bordery',       20,    INITOPT),
	    ('iconmargin',    20,    INITOPT),
	    ('iconpos',       None,  INITOPT),
	)
	self.defineoptions(kw, optiondefs)

	# Initialise the base class (after defining the options).
	Dialog.__init__(self, parent)

	# Create the components.
	interior = self.interior()

	self._message = self.createcomponent('message',
		(), None,
		Tkinter.Label, (interior,))

        iconpos = self['iconpos']
	iconmargin = self['iconmargin']
	borderx = self['borderx']
	bordery = self['bordery']
	border_right = 2
	border_bottom = 2
	if iconpos is None:
	    self._message.grid(column = 1, row = 1)
	else:
	    self._icon = self.createcomponent('icon',
		    (), None,
		    Tkinter.Label, (interior,))
	    if iconpos not in 'nsew':
		raise ValueError('bad iconpos option "%s":  should be n, s, e, or w'
                                 % iconpos)

	    if iconpos in 'nw':
		icon = 1
		message = 3
	    else:
		icon = 3
		message = 1

	    if iconpos in 'ns':
		# vertical layout
		self._icon.grid(column = 1, row = icon)
		self._message.grid(column = 1, row = message)
		interior.grid_rowconfigure(2, minsize = iconmargin)
		border_bottom = 4
	    else:
		# horizontal layout
		self._icon.grid(column = icon, row = 1)
		self._message.grid(column = message, row = 1)
		interior.grid_columnconfigure(2, minsize = iconmargin)
		border_right = 4

	interior.grid_columnconfigure(0, minsize = borderx)
	interior.grid_rowconfigure(0, minsize = bordery)
	interior.grid_columnconfigure(border_right, minsize = borderx)
	interior.grid_rowconfigure(border_bottom, minsize = bordery)


	# Check keywords and initialise options.
	self.initialiseoptions()

######################################################################
### File: PmwNoteBook.py
import string
import types
import Tkinter


class NoteBook(MegaArchetype):

    def __init__(self, parent = None, **kw):

	# Define the megawidget options.
	
        optiondefs = (
	    ('hull_highlightthickness',  0,           None),
	    ('hull_borderwidth',         0,           None),
            ('arrownavigation',          1,           INITOPT),
            ('borderwidth',              2,           INITOPT),
            ('createcommand',            None,        None),
            ('lowercommand',             None,        None),
            ('pagemargin',               4,           INITOPT),
            ('raisecommand',             None,        None),
	    ('tabpos',                   'n',         INITOPT),
        )
	self.defineoptions(kw, optiondefs, dynamicGroups = ('Page', 'Tab'))

	# Initialise the base class (after defining the options).
	MegaArchetype.__init__(self, parent, Tkinter.Canvas)

        self.bind('<Map>', self._handleMap)
        self.bind('<Configure>', self._handleConfigure)

        tabpos = self['tabpos']
	if tabpos is not None and tabpos != 'n':
            raise ValueError('bad tabpos option %s:  should be n or None'
                             % repr(tabpos))
        self._withTabs = (tabpos is not None)
        self._pageMargin = self['pagemargin']
        self._borderWidth = self['borderwidth']

        # Use a dictionary as a set of bits indicating what needs to
        # be redisplayed the next time _layout() is called.  If
        # dictionary contains 'topPage' key, the value is the new top
        # page to be displayed.  None indicates that all pages have
        # been deleted and that _layout() should draw a border under where
        # the tabs should be.
        self._pending = {}
        self._pending['size'] = 1
        self._pending['borderColor'] = 1
        self._pending['topPage'] = None
        if self._withTabs:
            self._pending['tabs'] = 1

        self._canvasSize = None       # This gets set by <Configure> events

        # Set initial height of space for tabs
        if self._withTabs:
            self.tabBottom = 35
        else:
            self.tabBottom = 0

        self._lightBorderColor, self._darkBorderColor = \
                Color.bordercolors(self, self['hull_background'])

	self._pageNames   = []        # List of page names

        # Map from page name to page info.  Each item is itself a
        # dictionary containing the following items:
        #   page           the Tkinter.Frame widget for the page
        #   created        set to true the first time the page is raised
        #   tabbutton      the Tkinter.Button widget for the button (if any)
        #   tabreqwidth    requested width of the tab
        #   tabreqheight   requested height of the tab
        #   tabitems       the canvas items for the button: the button
        #                  window item, the lightshadow and the darkshadow
        #   left           the left and right canvas coordinates of the tab
        #   right
	self._pageAttrs   = {}

        # Name of page currently on top (actually displayed, using
        # create_window, not pending).  Ignored if current top page
        # has been deleted or new top page is pending.  None indicates
        # no pages in notebook.
	self._topPageName = None

        # Canvas items used:
        #   Per tab: 
        #       top and left shadow
        #       right shadow
        #       button
        #   Per notebook: 
        #       page
        #       top page
        #       left shadow
        #       bottom and right shadow
        #       top (one or two items)

        # Canvas tags used:
        #   lighttag      - top and left shadows of tabs and page
        #   darktag       - bottom and right shadows of tabs and page
        #                   (if no tabs then these are reversed)
        #                   (used to color the borders by recolorborders)

        # Create page border shadows.
        if self._withTabs:
            self._pageLeftBorder = self.create_polygon(0, 0, 0, 0, 0, 0,
                fill = self._lightBorderColor, tags = 'lighttag')
            self._pageBottomRightBorder = self.create_polygon(0, 0, 0, 0, 0, 0,
                fill = self._darkBorderColor, tags = 'darktag')
            self._pageTop1Border = self.create_polygon(0, 0, 0, 0, 0, 0,
                fill = self._darkBorderColor, tags = 'lighttag')
            self._pageTop2Border = self.create_polygon(0, 0, 0, 0, 0, 0,
                fill = self._darkBorderColor, tags = 'lighttag')
        else:
            self._pageLeftBorder = self.create_polygon(0, 0, 0, 0, 0, 0,
                fill = self._darkBorderColor, tags = 'darktag')
            self._pageBottomRightBorder = self.create_polygon(0, 0, 0, 0, 0, 0,
                fill = self._lightBorderColor, tags = 'lighttag')
            self._pageTopBorder = self.create_polygon(0, 0, 0, 0, 0, 0,
                fill = self._darkBorderColor, tags = 'darktag')

	# Check keywords and initialise options.
	self.initialiseoptions()

    def insert(self, pageName, before = 0, **kw):
	if self._pageAttrs.has_key(pageName):
	    msg = 'Page "%s" already exists.' % pageName
	    raise ValueError(msg)

        # Do this early to catch bad <before> spec before creating any items.
	beforeIndex = self.index(before, 1)

        pageOptions = {}
        if self._withTabs:
            # Default tab button options.
            tabOptions = {
                'text' : pageName,
                'borderwidth' : 0,
            }

        # Divide the keyword options into the 'page_' and 'tab_' options.
        for key in kw.keys():
            if key[:5] == 'page_':
                pageOptions[key[5:]] = kw[key]
                del kw[key]
            elif self._withTabs and key[:4] == 'tab_':
                tabOptions[key[4:]] = kw[key]
                del kw[key]
            else:
		raise KeyError('Unknown option "' + key + '"')

        # Create the frame to contain the page.
	page = apply(self.createcomponent, (pageName,
		(), 'Page',
		Tkinter.Frame, self._hull), pageOptions)

        attributes = {}
        attributes['page'] = page
        attributes['created'] = 0

        if self._withTabs:
            # Create the button for the tab.
            def raiseThisPage(self = self, pageName = pageName):
                self.selectpage(pageName)
            tabOptions['command'] = raiseThisPage
            tab = apply(self.createcomponent, (pageName + '-tab',
                    (), 'Tab',
                    Tkinter.Button, self._hull), tabOptions)

            if self['arrownavigation']:
                # Allow the use of the arrow keys for Tab navigation:
                def next(event, self = self, pageName = pageName):
                    self.nextpage(pageName)
                def prev(event, self = self, pageName = pageName):
                    self.previouspage(pageName)
                tab.bind('<Left>', prev)
                tab.bind('<Right>', next)

            attributes['tabbutton'] = tab
            attributes['tabreqwidth'] = tab.winfo_reqwidth()
            attributes['tabreqheight'] = tab.winfo_reqheight()

            # Create the canvas item to manage the tab's button and the items
            # for the tab's shadow.
            windowitem = self.create_window(0, 0, window = tab, anchor = 'nw')
            lightshadow = self.create_polygon(0, 0, 0, 0, 0, 0,
                tags = 'lighttag', fill = self._lightBorderColor)
            darkshadow = self.create_polygon(0, 0, 0, 0, 0, 0,
                tags = 'darktag', fill = self._darkBorderColor)
            attributes['tabitems'] = (windowitem, lightshadow, darkshadow)
            self._pending['tabs'] = 1

        self._pageAttrs[pageName] = attributes
	self._pageNames.insert(beforeIndex, pageName)

        # If this is the first page added, make it the new top page
        # and call the create and raise callbacks.
        if self.getcurselection() is None:
            self._pending['topPage'] = pageName
            self._raiseNewTop(pageName)

        self._layout()
        return page
  		
    def add(self, pageName, **kw):
        return apply(self.insert, (pageName, len(self._pageNames)), kw)

    def delete(self, *pageNames):
        newTopPage = 0
        for page in pageNames:
            pageIndex = self.index(page)
            pageName = self._pageNames[pageIndex]
            pageInfo = self._pageAttrs[pageName]

            if self.getcurselection() == pageName:
                if len(self._pageNames) == 1:
                    newTopPage = 0
                    self._pending['topPage'] = None
                elif pageIndex == len(self._pageNames) - 1:
                    newTopPage = 1
                    self._pending['topPage'] = self._pageNames[pageIndex - 1]
                else:
                    newTopPage = 1
                    self._pending['topPage'] = self._pageNames[pageIndex + 1]

            if self._topPageName == pageName:
                self._hull.delete(self._topPageItem)
                self._topPageName = None
                                
            if self._withTabs:
                self.destroycomponent(pageName + '-tab')
                apply(self._hull.delete, pageInfo['tabitems'])
            self.destroycomponent(pageName)
            del self._pageAttrs[pageName]
            del self._pageNames[pageIndex]

        # If the old top page was deleted and there are still pages
        # left in the notebook, call the create and raise callbacks.
        if newTopPage:
            pageName = self._pending['topPage']
            self._raiseNewTop(pageName)

        if self._withTabs:
            self._pending['tabs'] = 1
        self._layout()

    def page(self, pageIndex):
        pageName = self._pageNames[self.index(pageIndex)]
	return self._pageAttrs[pageName]['page']

    def pagenames(self):
	return list(self._pageNames)

    def getcurselection(self):
        if self._pending.has_key('topPage'):
            return self._pending['topPage']
        else:
            return self._topPageName

    def tab(self, pageIndex):
        if self._withTabs:
            pageName = self._pageNames[self.index(pageIndex)]
            return self._pageAttrs[pageName]['tabbutton']
        else:
            return None

    def index(self, index, forInsert = 0):
	listLength = len(self._pageNames)
        if isinstance(index, int):
	    if forInsert and index <= listLength:
		return index
	    elif not forInsert and index < listLength:
		return index
	    else:
		raise ValueError('index "%s" is out of range' % index)
	elif index is END:
	    if forInsert:
		return listLength
	    elif listLength > 0:
		return listLength - 1
	    else:
		raise ValueError('NoteBook has no pages')
	elif index is SELECT:
	    if listLength == 0:
		raise ValueError('NoteBook has no pages')
            return self._pageNames.index(self.getcurselection())
	else:
            if index in self._pageNames:
                return self._pageNames.index(index)
	    validValues = 'a name, a number, END or SELECT'
	    raise ValueError('bad index "%s": must be %s' % (index, validValues))

    def selectpage(self, page):
        pageName = self._pageNames[self.index(page)]
        oldTopPage = self.getcurselection()
        if pageName != oldTopPage:
            self._pending['topPage'] = pageName
            if oldTopPage == self._topPageName:
                self._hull.delete(self._topPageItem)
            cmd = self['lowercommand']
            if cmd is not None:
                cmd(oldTopPage)
            self._raiseNewTop(pageName)

            self._layout()

        # Set focus to the tab of new top page:
        if self._withTabs and self['arrownavigation']:
            self._pageAttrs[pageName]['tabbutton'].focus_set()

    def previouspage(self, pageIndex = None):
        if pageIndex is None:
            curpage = self.index(SELECT)
        else:
            curpage = self.index(pageIndex)
	if curpage > 0:
	    self.selectpage(curpage - 1)

    def nextpage(self, pageIndex = None):
        if pageIndex is None:
            curpage = self.index(SELECT)
        else:
            curpage = self.index(pageIndex)
	if curpage < len(self._pageNames) - 1:
	    self.selectpage(curpage + 1)

    def setnaturalsize(self, pageNames = None):
        self.update_idletasks()
        maxPageWidth = 1
        maxPageHeight = 1
        if pageNames is None:
            pageNames = self.pagenames()
        for pageName in pageNames:
            pageInfo = self._pageAttrs[pageName]
            page = pageInfo['page']
            w = page.winfo_reqwidth()
            h = page.winfo_reqheight()
            if maxPageWidth < w:
                maxPageWidth = w
            if maxPageHeight < h:
                maxPageHeight = h
        pageBorder = self._borderWidth + self._pageMargin
        width = maxPageWidth + pageBorder * 2
        height = maxPageHeight + pageBorder * 2

        if self._withTabs:
            maxTabHeight = 0
            for pageInfo in self._pageAttrs.values():
                if maxTabHeight < pageInfo['tabreqheight']:
                    maxTabHeight = pageInfo['tabreqheight']
            height = height + maxTabHeight + self._borderWidth * 1.5

        # Note that, since the hull is a canvas, the width and height
        # options specify the geometry *inside* the borderwidth and
        # highlightthickness.
        self.configure(hull_width = width, hull_height = height)

    def recolorborders(self):
        self._pending['borderColor'] = 1
        self._layout()

    def _handleMap(self, event):
        self._layout()

    def _handleConfigure(self, event):
        self._canvasSize = (event.width, event.height)
        self._pending['size'] = 1
        self._layout()

    def _raiseNewTop(self, pageName):
        if not self._pageAttrs[pageName]['created']:
            self._pageAttrs[pageName]['created'] = 1
            cmd = self['createcommand']
            if cmd is not None:
                cmd(pageName)
        cmd = self['raisecommand']
        if cmd is not None:
            cmd(pageName)

    # This is the vertical layout of the notebook, from top (assuming
    # tabpos is 'n'):
    #     hull highlightthickness (top)
    #     hull borderwidth (top)
    #     borderwidth (top border of tabs)
    #     borderwidth * 0.5 (space for bevel)
    #     tab button (maximum of requested height of all tab buttons)
    #     borderwidth (border between tabs and page)
    #     pagemargin (top)
    #     the page itself
    #     pagemargin (bottom)
    #     borderwidth (border below page)
    #     hull borderwidth (bottom)
    #     hull highlightthickness (bottom)
    #
    # canvasBorder is sum of top two elements.
    # tabBottom is sum of top five elements.
    #
    # Horizontal layout (and also vertical layout when tabpos is None):
    #     hull highlightthickness
    #     hull borderwidth
    #     borderwidth
    #     pagemargin
    #     the page itself
    #     pagemargin
    #     borderwidth
    #     hull borderwidth
    #     hull highlightthickness
    #
    def _layout(self):
        if not self.winfo_ismapped() or self._canvasSize is None:
            # Don't layout if the window is not displayed, or we
            # haven't yet received a <Configure> event.
            return

        hullWidth, hullHeight = self._canvasSize
        borderWidth = self._borderWidth
        canvasBorder = string.atoi(self._hull['borderwidth']) + \
            string.atoi(self._hull['highlightthickness'])
        if not self._withTabs:
            self.tabBottom = canvasBorder
        oldTabBottom = self.tabBottom

        if self._pending.has_key('borderColor'):
            self._lightBorderColor, self._darkBorderColor = \
                    Color.bordercolors(self, self['hull_background'])

        # Draw all the tabs.
        if self._withTabs and (self._pending.has_key('tabs') or
                self._pending.has_key('size')):
            # Find total requested width and maximum requested height
            # of tabs.
            sumTabReqWidth = 0
            maxTabHeight = 0
            for pageInfo in self._pageAttrs.values():
                sumTabReqWidth = sumTabReqWidth + pageInfo['tabreqwidth']
                if maxTabHeight < pageInfo['tabreqheight']:
                    maxTabHeight = pageInfo['tabreqheight']
            if maxTabHeight != 0:
                # Add the top tab border plus a bit for the angled corners
                self.tabBottom = canvasBorder + maxTabHeight + borderWidth * 1.5

            # Prepare for drawing the border around each tab button.
            tabTop = canvasBorder
            tabTop2 = tabTop + borderWidth
            tabTop3 = tabTop + borderWidth * 1.5
            tabBottom2 = self.tabBottom
            tabBottom = self.tabBottom + borderWidth

            numTabs = len(self._pageNames)
            availableWidth = hullWidth - 2 * canvasBorder - \
                numTabs * 2 * borderWidth
            x = canvasBorder
            cumTabReqWidth = 0
            cumTabWidth = 0

            # Position all the tabs.
            for pageName in self._pageNames:
                pageInfo = self._pageAttrs[pageName]
                (windowitem, lightshadow, darkshadow) = pageInfo['tabitems']
                if sumTabReqWidth <= availableWidth:
                    tabwidth = pageInfo['tabreqwidth']
                else:
                    # This ugly calculation ensures that, when the
                    # notebook is not wide enough for the requested
                    # widths of the tabs, the total width given to
                    # the tabs exactly equals the available width,
                    # without rounding errors.
                    cumTabReqWidth = cumTabReqWidth + pageInfo['tabreqwidth']
                    tmp = (2*cumTabReqWidth*availableWidth + sumTabReqWidth) \
                            / (2 * sumTabReqWidth)
                    tabwidth = tmp - cumTabWidth
                    cumTabWidth = tmp

                # Position the tab's button canvas item.
                self.coords(windowitem, x + borderWidth, tabTop3)
                self.itemconfigure(windowitem,
                    width = tabwidth, height = maxTabHeight)

                # Make a beautiful border around the tab.
                left = x
                left2 = left + borderWidth
                left3 = left + borderWidth * 1.5
                right = left + tabwidth + 2 * borderWidth
                right2 = left + tabwidth + borderWidth
                right3 = left + tabwidth + borderWidth * 0.5

                self.coords(lightshadow, 
                    left, tabBottom2, left, tabTop2, left2, tabTop,
                    right2, tabTop, right3, tabTop2, left3, tabTop2,
                    left2, tabTop3, left2, tabBottom,
                    )
                self.coords(darkshadow, 
                    right2, tabTop, right, tabTop2, right, tabBottom2,
                    right2, tabBottom, right2, tabTop3, right3, tabTop2,
                    )
                pageInfo['left'] = left
                pageInfo['right'] = right

                x = x + tabwidth + 2 * borderWidth

        # Redraw shadow under tabs so that it appears that tab for old
        # top page is lowered and that tab for new top page is raised.
        if self._withTabs and (self._pending.has_key('topPage') or
                self._pending.has_key('tabs') or self._pending.has_key('size')):

            if self.getcurselection() is None:
                # No pages, so draw line across top of page area.
                self.coords(self._pageTop1Border,
                    canvasBorder, self.tabBottom,
                    hullWidth - canvasBorder, self.tabBottom,
                    hullWidth - canvasBorder - borderWidth,
                        self.tabBottom + borderWidth,
                    borderWidth + canvasBorder, self.tabBottom + borderWidth,
                    )

                # Ignore second top border.
                self.coords(self._pageTop2Border, 0, 0, 0, 0, 0, 0)
            else:
                # Draw two lines, one on each side of the tab for the
                # top page, so that the tab appears to be raised.
                pageInfo = self._pageAttrs[self.getcurselection()]
                left = pageInfo['left']
                right = pageInfo['right']
                self.coords(self._pageTop1Border,
                    canvasBorder, self.tabBottom,
                    left, self.tabBottom,
                    left + borderWidth, self.tabBottom + borderWidth,
                    canvasBorder + borderWidth, self.tabBottom + borderWidth,
                    )

                self.coords(self._pageTop2Border,
                    right, self.tabBottom,
                    hullWidth - canvasBorder, self.tabBottom,
                    hullWidth - canvasBorder - borderWidth,
                        self.tabBottom + borderWidth,
                    right - borderWidth, self.tabBottom + borderWidth,
                    )

            # Prevent bottom of dark border of tabs appearing over
            # page top border.
            self.tag_raise(self._pageTop1Border)
            self.tag_raise(self._pageTop2Border)

        # Position the page border shadows.
        if self._pending.has_key('size') or oldTabBottom != self.tabBottom:

            self.coords(self._pageLeftBorder,
                canvasBorder, self.tabBottom,
                borderWidth + canvasBorder,
                    self.tabBottom + borderWidth,
                borderWidth + canvasBorder,
                    hullHeight - canvasBorder - borderWidth,
                canvasBorder, hullHeight - canvasBorder,
                )

            self.coords(self._pageBottomRightBorder,
                hullWidth - canvasBorder, self.tabBottom,
                hullWidth - canvasBorder, hullHeight - canvasBorder,
                canvasBorder, hullHeight - canvasBorder,
                borderWidth + canvasBorder,
                    hullHeight - canvasBorder - borderWidth,
                hullWidth - canvasBorder - borderWidth,
                    hullHeight - canvasBorder - borderWidth,
                hullWidth - canvasBorder - borderWidth,
                    self.tabBottom + borderWidth,
                )

            if not self._withTabs:
                self.coords(self._pageTopBorder,
                    canvasBorder, self.tabBottom,
                    hullWidth - canvasBorder, self.tabBottom,
                    hullWidth - canvasBorder - borderWidth,
                        self.tabBottom + borderWidth,
                    borderWidth + canvasBorder, self.tabBottom + borderWidth,
                    )

        # Color borders.
        if self._pending.has_key('borderColor'):
            self.itemconfigure('lighttag', fill = self._lightBorderColor)
            self.itemconfigure('darktag', fill = self._darkBorderColor)

        newTopPage = self._pending.get('topPage')
        pageBorder = borderWidth + self._pageMargin

        # Raise new top page.
        if newTopPage is not None:
            self._topPageName = newTopPage
            self._topPageItem = self.create_window(
                pageBorder + canvasBorder, self.tabBottom + pageBorder,
                window = self._pageAttrs[newTopPage]['page'],
                anchor = 'nw',
                )

        # Change position of top page if tab height has changed.
        if self._topPageName is not None and oldTabBottom != self.tabBottom:
            self.coords(self._topPageItem,
                    pageBorder + canvasBorder, self.tabBottom + pageBorder)

        # Change size of top page if,
        #   1) there is a new top page.
        #   2) canvas size has changed, but not if there is no top
        #      page (eg:  initially or when all pages deleted).
        #   3) tab height has changed, due to difference in the height of a tab
        if (newTopPage is not None or \
                self._pending.has_key('size') and self._topPageName is not None
                or oldTabBottom != self.tabBottom):
            self.itemconfigure(self._topPageItem,
                width = hullWidth - 2 * canvasBorder - pageBorder * 2,
                height = hullHeight - 2 * canvasBorder - pageBorder * 2 -
                    (self.tabBottom - canvasBorder),
                )

        self._pending = {}

# Need to do forwarding to get the pack, grid, etc methods. 
# Unfortunately this means that all the other canvas methods are also
# forwarded.
forwardmethods(NoteBook, Tkinter.Canvas, '_hull')

######################################################################
### File: PmwOptionMenu.py
import types
import Tkinter


class OptionMenu(MegaWidget):

    def __init__(self, parent = None, **kw):

	# Define the megawidget options.
	
	optiondefs = (
	    ('command',        None,       None),
            ('items',          (),         INITOPT),
            ('initialitem',    None,       INITOPT),
	    ('labelmargin',    0,          INITOPT),
	    ('labelpos',       None,       INITOPT),
	    ('sticky',         'ew',       INITOPT),
	)
	self.defineoptions(kw, optiondefs)

	# Initialise the base class (after defining the options).
	MegaWidget.__init__(self, parent)

	# Create the components.
	interior = self.interior()

	self._menubutton = self.createcomponent('menubutton',
		(), None,
		Tkinter.Menubutton, (interior,),
		borderwidth = 2,
		indicatoron = 1,
		relief = 'raised',
		anchor = 'c',
		highlightthickness = 2,
		direction = 'flush',
                takefocus = 1,
	)
	self._menubutton.grid(column = 2, row = 2, sticky = self['sticky'])

	self._menu = self.createcomponent('menu',
		(), None,
		Tkinter.Menu, (self._menubutton,),
		tearoff=0
	)
	self._menubutton.configure(menu = self._menu)

	interior.grid_columnconfigure(2, weight = 1)
	interior.grid_rowconfigure(2, weight = 1)

        # Create the label.
        self.createlabel(interior)

        # Add the items specified by the initialisation option.
	self._itemList = []
        self.setitems(self['items'], self['initialitem'])

	# Check keywords and initialise options.
	self.initialiseoptions()

    def setitems(self, items, index = None):

        # Clean up old items and callback commands.
        for oldIndex in range(len(self._itemList)):
            tclCommandName = str(self._menu.entrycget(oldIndex, 'command'))
            if tclCommandName != '':   
                self._menu.deletecommand(tclCommandName)
        self._menu.delete(0, 'end')
	self._itemList = list(items)

	# Set the items in the menu component.
        for item in items:
            self._menu.add_command(label = item,
		command = lambda self = self, item = item: self._invoke(item))

	# Set the currently selected value.
	if index is None:
            var = str(self._menubutton.cget('textvariable'))
	    if var != '':
		# None means do not change text variable.
		return
	    if len(items) == 0:
		text = ''
	    elif str(self._menubutton.cget('text')) in items:
                # Do not change selection if it is still valid
		return
	    else:
		text = items[0]
	else:
	    index = self.index(index)
	    text = self._itemList[index]

        self.setvalue(text)

    def getcurselection(self):
	var = str(self._menubutton.cget('textvariable'))
	if var == '':
	    return str(self._menubutton.cget('text'))
	else:
	    return self._menu.tk.globalgetvar(var)

    def getvalue(self):
        return self.getcurselection()

    def setvalue(self, text):
	var = str(self._menubutton.cget('textvariable'))
	if var == '':
	    self._menubutton.configure(text = text)
	else:
	    self._menu.tk.globalsetvar(var, text)

    def index(self, index):
	listLength = len(self._itemList)
        if isinstance(index, int):
	    if index < listLength:
		return index
	    else:
		raise ValueError('index "%s" is out of range' % index)
	elif index is END:
	    if listLength > 0:
		return listLength - 1
	    else:
		raise ValueError('OptionMenu has no items')
	else:
	    if index is SELECT:
		if listLength > 0:
		    index = self.getcurselection()
		else:
		    raise ValueError('OptionMenu has no items')
            if index in self._itemList:
                return self._itemList.index(index)
	    raise ValueError('bad index "%s": must be a name, a number, END or SELECT'
                             % (index,))

    def invoke(self, index = SELECT):
	index = self.index(index)
	text = self._itemList[index]

        return self._invoke(text)

    def _invoke(self, text):
        self.setvalue(text)

	command = self['command']
	if callable(command):
	    return command(text)

######################################################################
### File: PmwPanedWidget.py
# PanedWidget
# a frame which may contain several resizable sub-frames

import string
import sys
import types
import Tkinter


class PanedWidget(MegaWidget):

    def __init__(self, parent = None, **kw):

	# Define the megawidget options.
	
	optiondefs = (
            ('command',            None,         None),
            ('orient',             'vertical',   INITOPT),
            ('separatorrelief',    'sunken',     INITOPT),
            ('separatorthickness', 2,            INITOPT),
            ('handlesize',         8,            INITOPT),
            ('hull_width',         400,          None),
            ('hull_height',        400,          None),
	)
	self.defineoptions(kw, optiondefs,
                dynamicGroups = ('Frame', 'Separator', 'Handle'))

	# Initialise the base class (after defining the options).
	MegaWidget.__init__(self, parent)

	self.bind('<Configure>', self._handleConfigure)

	if self['orient'] not in ('horizontal', 'vertical'):
	    raise ValueError("bad orient option %s: must be either 'horizontal' or "
                             "'vertical'" % repr(self['orient']))

        self._separatorThickness = self['separatorthickness']
        self._handleSize = self['handlesize']
	self._paneNames = []            # List of pane names
	self._paneAttrs = {}            # Map from pane name to pane info

	self._timerId = None
	self._frame = {}
	self._separator = []
	self._button = []
	self._totalSize = 0
	self._movePending = 0
	self._relsize = {}
	self._relmin = {}
	self._relmax = {}
	self._size = {}
	self._min = {}
	self._max = {}
	self._rootp = None
	self._curSize = None
	self._beforeLimit = None
	self._afterLimit = None
	self._buttonIsDown = 0
	self._majorSize = 100
	self._minorSize = 100

	# Check keywords and initialise options.
	self.initialiseoptions()

    def insert(self, name, before = 0, **kw):
	# Parse <kw> for options.
        self._initPaneOptions(name)
	self._parsePaneOptions(name, kw)

	insertPos = self._nameToIndex(before)
	atEnd = (insertPos == len(self._paneNames))

	# Add the frame.
	self._paneNames[insertPos:insertPos] = [name]
	self._frame[name] = self.createcomponent(name,
		(), 'Frame',
		Tkinter.Frame, (self.interior(),))

	# Add separator, if necessary.
	if len(self._paneNames) > 1:
	    self._addSeparator()
	else:
	    self._separator.append(None)
	    self._button.append(None)

	# Add the new frame and adjust the PanedWidget
	if atEnd:
	    size = self._size[name]
	    if size > 0 or self._relsize[name] is not None:
		if self['orient'] == 'vertical':
		    self._frame[name].place(x=0, relwidth=1,
					    height=size, y=self._totalSize)
		else:
		    self._frame[name].place(y=0, relheight=1,
					    width=size, x=self._totalSize)
	    else:
		if self['orient'] == 'vertical':
		    self._frame[name].place(x=0, relwidth=1,
					    y=self._totalSize)
		else:
		    self._frame[name].place(y=0, relheight=1,
					    x=self._totalSize)
	else:
	    self._updateSizes()

	self._totalSize = self._totalSize + self._size[name]
	return self._frame[name]

    def add(self, name, **kw):
        return apply(self.insert, (name, len(self._paneNames)), kw)

    def delete(self, name):
	deletePos = self._nameToIndex(name)
	name = self._paneNames[deletePos]
	self.destroycomponent(name)
	del self._paneNames[deletePos]
	del self._frame[name]
	del self._size[name]
	del self._min[name]
	del self._max[name]
	del self._relsize[name]
	del self._relmin[name]
	del self._relmax[name]

	last = len(self._paneNames)
	del self._separator[last]
	del self._button[last]
        if last > 0:
            self.destroycomponent(self._sepName(last))
            self.destroycomponent(self._buttonName(last))

	self._plotHandles()

    def setnaturalsize(self):
        self.update_idletasks()
        totalWidth = 0
        totalHeight = 0
        maxWidth = 0
        maxHeight = 0
	for name in self._paneNames:
            frame = self._frame[name]
            w = frame.winfo_reqwidth()
            h = frame.winfo_reqheight()
            totalWidth = totalWidth + w
            totalHeight = totalHeight + h
            if maxWidth < w:
                maxWidth = w
            if maxHeight < h:
                maxHeight = h

        # Note that, since the hull is a frame, the width and height
        # options specify the geometry *outside* the borderwidth and
        # highlightthickness.
        bw = string.atoi(str(self.cget('hull_borderwidth')))
        hl = string.atoi(str(self.cget('hull_highlightthickness')))
        extra = (bw + hl) * 2
        if str(self.cget('orient')) == 'horizontal':
            totalWidth = totalWidth + extra
            maxHeight = maxHeight + extra
            self.configure(hull_width = totalWidth, hull_height = maxHeight)
        else:
            totalHeight = (totalHeight + extra +
                    (len(self._paneNames) - 1) * self._separatorThickness)
            maxWidth = maxWidth + extra
            self.configure(hull_width = maxWidth, hull_height = totalHeight)

    def move(self, name, newPos, newPosOffset = 0):

        # see if we can spare ourselves some work
        numPanes = len(self._paneNames)
        if numPanes < 2:
            return

        newPos = self._nameToIndex(newPos) + newPosOffset
        if newPos < 0 or newPos >=numPanes:
            return

        deletePos = self._nameToIndex(name)

        if deletePos == newPos:
            # inserting over ourself is a no-op
            return

        # delete name from old position in list
        name = self._paneNames[deletePos]
        del self._paneNames[deletePos]

        # place in new position
        self._paneNames[newPos:newPos] = [name]

        # force everything to redraw
        self._plotHandles()
        self._updateSizes()

    def _nameToIndex(self, nameOrIndex):
	try:
	    pos = self._paneNames.index(nameOrIndex)
	except ValueError:
	    pos = nameOrIndex

	return pos

    def _initPaneOptions(self, name):
	# Set defaults.
	self._size[name] = 0
	self._relsize[name] = None
	self._min[name] = 0
	self._relmin[name] = None
	self._max[name] = 100000
	self._relmax[name] = None

    def _parsePaneOptions(self, name, args):
	# Parse <args> for options.
	for arg, value in args.items():
            if isinstance(value, float):
		relvalue = value
		value = self._absSize(relvalue)
	    else:
		relvalue = None

	    if arg == 'size':
		self._size[name], self._relsize[name] = value, relvalue
	    elif arg == 'min':
		self._min[name], self._relmin[name] = value, relvalue
	    elif arg == 'max':
		self._max[name], self._relmax[name] = value, relvalue
	    else:
		raise ValueError('keyword must be "size", "min", or "max"')

    def _absSize(self, relvalue):
	return int(round(relvalue * self._majorSize))

    def _sepName(self, n):
	return 'separator-%d' % n

    def _buttonName(self, n):
	return 'handle-%d' % n

    def _addSeparator(self):
	n = len(self._paneNames) - 1

	downFunc = lambda event, s = self, num=n: s._btnDown(event, num)
	upFunc = lambda event, s = self, num=n: s._btnUp(event, num)
	moveFunc = lambda event, s = self, num=n: s._btnMove(event, num)

	# Create the line dividing the panes.
	sep = self.createcomponent(self._sepName(n),
		(), 'Separator',
		Tkinter.Frame, (self.interior(),),
		borderwidth = 1,
		relief = self['separatorrelief'])
	self._separator.append(sep)

	sep.bind('<ButtonPress-1>', downFunc)
	sep.bind('<Any-ButtonRelease-1>', upFunc)
	sep.bind('<B1-Motion>', moveFunc)

	if self['orient'] == 'vertical':
	    cursor = 'sb_v_double_arrow'
	    sep.configure(height = self._separatorThickness,
                    width = 10000, cursor = cursor)
	else:
	    cursor = 'sb_h_double_arrow'
	    sep.configure(width = self._separatorThickness,
                    height = 10000, cursor = cursor)

	self._totalSize = self._totalSize + self._separatorThickness

	# Create the handle on the dividing line.
	handle = self.createcomponent(self._buttonName(n),
		(), 'Handle',
		Tkinter.Frame, (self.interior(),),
		    relief = 'raised',
		    borderwidth = 1,
		    width = self._handleSize,
		    height = self._handleSize,
		    cursor = cursor,
		)
	self._button.append(handle)

	handle.bind('<ButtonPress-1>', downFunc)
	handle.bind('<Any-ButtonRelease-1>', upFunc)
	handle.bind('<B1-Motion>', moveFunc)

	self._plotHandles()

	for i in range(1, len(self._paneNames)):
	    self._separator[i].tkraise()
	for i in range(1, len(self._paneNames)):
	    self._button[i].tkraise()

    def _btnUp(self, event, item):
	self._buttonIsDown = 0
	self._updateSizes()
	try:
	    self._button[item].configure(relief='raised')
	except:
	    pass

    def _btnDown(self, event, item):
	self._button[item].configure(relief='sunken')
	self._getMotionLimit(item)
	self._buttonIsDown = 1
	self._movePending = 0

    def _handleConfigure(self, event = None):
	self._getNaturalSizes()
	if self._totalSize == 0:
	    return

	iterRange = list(self._paneNames)
	iterRange.reverse()
	if self._majorSize > self._totalSize:
	    n = self._majorSize - self._totalSize
	    self._iterate(iterRange, self._grow, n)
	elif self._majorSize < self._totalSize:
	    n = self._totalSize - self._majorSize
	    self._iterate(iterRange, self._shrink, n)

	self._plotHandles()
	self._updateSizes()

    def _getNaturalSizes(self):
	# Must call this in order to get correct winfo_width, winfo_height
	self.update_idletasks()

	self._totalSize = 0

	if self['orient'] == 'vertical':
	    self._majorSize = self.winfo_height()
	    self._minorSize = self.winfo_width()
	    majorspec = Tkinter.Frame.winfo_reqheight
	else:
	    self._majorSize = self.winfo_width()
	    self._minorSize = self.winfo_height()
	    majorspec = Tkinter.Frame.winfo_reqwidth

        bw = string.atoi(str(self.cget('hull_borderwidth')))
        hl = string.atoi(str(self.cget('hull_highlightthickness')))
        extra = (bw + hl) * 2
        self._majorSize = self._majorSize - extra
        self._minorSize = self._minorSize - extra

	if self._majorSize < 0:
	    self._majorSize = 0
	if self._minorSize < 0:
	    self._minorSize = 0

	for name in self._paneNames:
	    # adjust the absolute sizes first...
	    if self._relsize[name] is None:
		#special case
		if self._size[name] == 0:
		    self._size[name] = apply(majorspec, (self._frame[name],))
		    self._setrel(name)
	    else:
		self._size[name] = self._absSize(self._relsize[name])

	    if self._relmin[name] is not None:
		self._min[name] = self._absSize(self._relmin[name])
	    if self._relmax[name] is not None:
		self._max[name] = self._absSize(self._relmax[name])

	    # now adjust sizes
	    if self._size[name] < self._min[name]:
		self._size[name] = self._min[name]
		self._setrel(name)

	    if self._size[name] > self._max[name]:
		self._size[name] = self._max[name]
		self._setrel(name)

	    self._totalSize = self._totalSize + self._size[name]

	# adjust for separators
	self._totalSize = (self._totalSize +
                (len(self._paneNames) - 1) * self._separatorThickness)

    def _setrel(self, name):
	if self._relsize[name] is not None:
	    if self._majorSize != 0:
		self._relsize[name] = round(self._size[name]) / self._majorSize

    def _iterate(self, names, proc, n):
	for i in names:
	    n = apply(proc, (i, n))
	    if n == 0:
		break

    def _grow(self, name, n):
	canGrow = self._max[name] - self._size[name]

	if canGrow > n:
	    self._size[name] = self._size[name] + n
	    self._setrel(name)
	    return 0
	elif canGrow > 0:
	    self._size[name] = self._max[name]
	    self._setrel(name)
	    n = n - canGrow

	return n

    def _shrink(self, name, n):
	canShrink = self._size[name] - self._min[name]

	if canShrink > n:
	    self._size[name] = self._size[name] - n
	    self._setrel(name)
	    return 0
	elif canShrink > 0:
	    self._size[name] = self._min[name]
	    self._setrel(name)
	    n = n - canShrink

	return n

    def _updateSizes(self):
	totalSize = 0

	for name in self._paneNames:
	    size = self._size[name]
	    if self['orient'] == 'vertical':
		self._frame[name].place(x = 0, relwidth = 1,
					y = totalSize,
					height = size)
	    else:
		self._frame[name].place(y = 0, relheight = 1,
					x = totalSize,
					width = size)

	    totalSize = totalSize + size + self._separatorThickness

	# Invoke the callback command
	cmd = self['command']
	if callable(cmd):
	    cmd(map(lambda x, s = self: s._size[x], self._paneNames))

    def _plotHandles(self):
	if len(self._paneNames) == 0:
	    return

	if self['orient'] == 'vertical':
	    btnp = self._minorSize - 13
	else:
	    h = self._minorSize

	    if h > 18:
		btnp = 9
	    else:
		btnp = h - 9

	firstPane = self._paneNames[0]
	totalSize = self._size[firstPane]

	first = 1
	last = len(self._paneNames) - 1

	# loop from first to last, inclusive
	for i in range(1, last + 1):

	    handlepos = totalSize - 3
	    prevSize = self._size[self._paneNames[i - 1]]
	    nextSize = self._size[self._paneNames[i]]

	    offset1 = 0

	    if i == first:
		if prevSize < 4:
		    offset1 = 4 - prevSize
	    else:
		if prevSize < 8:
		    offset1 = (8 - prevSize) / 2

	    offset2 = 0

	    if i == last:
		if nextSize < 4:
		    offset2 = nextSize - 4
	    else:
		if nextSize < 8:
		    offset2 = (nextSize - 8) / 2

	    handlepos = handlepos + offset1

	    if self['orient'] == 'vertical':
		height = 8 - offset1 + offset2

		if height > 1:
		    self._button[i].configure(height = height)
		    self._button[i].place(x = btnp, y = handlepos)
		else:
		    self._button[i].place_forget()

		self._separator[i].place(x = 0, y = totalSize,
					 relwidth = 1)
	    else:
		width = 8 - offset1 + offset2

		if width > 1:
		    self._button[i].configure(width = width)
		    self._button[i].place(y = btnp, x = handlepos)
		else:
		    self._button[i].place_forget()

		self._separator[i].place(y = 0, x = totalSize,
					 relheight = 1)

	    totalSize = totalSize + nextSize + self._separatorThickness

    def pane(self, name):
	return self._frame[self._paneNames[self._nameToIndex(name)]]

    # Return the name of all panes
    def panes(self):
	return list(self._paneNames)

    def configurepane(self, name, **kw):
	name = self._paneNames[self._nameToIndex(name)]
	self._parsePaneOptions(name, kw)
	self._handleConfigure()

    def updatelayout(self):
	self._handleConfigure()

    def _getMotionLimit(self, item):
	curBefore = (item - 1) * self._separatorThickness
	minBefore, maxBefore = curBefore, curBefore

	for name in self._paneNames[:item]:
	    curBefore = curBefore + self._size[name]
	    minBefore = minBefore + self._min[name]
	    maxBefore = maxBefore + self._max[name]

	curAfter = (len(self._paneNames) - item) * self._separatorThickness
	minAfter, maxAfter = curAfter, curAfter
	for name in self._paneNames[item:]:
	    curAfter = curAfter + self._size[name]
	    minAfter = minAfter + self._min[name]
	    maxAfter = maxAfter + self._max[name]

	beforeToGo = min(curBefore - minBefore, maxAfter - curAfter)
	afterToGo = min(curAfter - minAfter, maxBefore - curBefore)

	self._beforeLimit = curBefore - beforeToGo
	self._afterLimit = curBefore + afterToGo
	self._curSize = curBefore

	self._plotHandles()

    # Compress the motion so that update is quick even on slow machines
    #
    # theRootp = root position (either rootx or rooty)
    def _btnMove(self, event, item):
	self._rootp = event

	if self._movePending == 0:
	    self._timerId = self.after_idle(
		    lambda s = self, i = item: s._btnMoveCompressed(i))
	    self._movePending = 1

    def destroy(self):
        if self._timerId is not None:
          self.after_cancel(self._timerId)
	  self._timerId = None
        MegaWidget.destroy(self)

    def _btnMoveCompressed(self, item):
	if not self._buttonIsDown:
	    return

	if self['orient'] == 'vertical':
	    p = self._rootp.y_root - self.winfo_rooty()
	else:
	    p = self._rootp.x_root - self.winfo_rootx()

	if p == self._curSize:
	    self._movePending = 0
	    return

	if p < self._beforeLimit:
	    p = self._beforeLimit

	if p >= self._afterLimit:
	    p = self._afterLimit

	self._calculateChange(item, p)
	self.update_idletasks()
	self._movePending = 0

    # Calculate the change in response to mouse motions
    def _calculateChange(self, item, p):

	if p < self._curSize:
	    self._moveBefore(item, p)
	elif p > self._curSize:
	    self._moveAfter(item, p)

	self._plotHandles()

    def _moveBefore(self, item, p):
	n = self._curSize - p

	# Shrink the frames before
	iterRange = list(self._paneNames[:item])
	iterRange.reverse()
	self._iterate(iterRange, self._shrink, n)

	# Adjust the frames after
	iterRange = self._paneNames[item:]
	self._iterate(iterRange, self._grow, n)

	self._curSize = p

    def _moveAfter(self, item, p):
	n = p - self._curSize

	# Shrink the frames after
	iterRange = self._paneNames[item:]
	self._iterate(iterRange, self._shrink, n)

	# Adjust the frames before
	iterRange = list(self._paneNames[:item])
	iterRange.reverse()
	self._iterate(iterRange, self._grow, n)

	self._curSize = p

######################################################################
### File: PmwPromptDialog.py
# Based on iwidgets2.2.0/promptdialog.itk code.



# A Dialog with an entryfield

class PromptDialog(Dialog):
    def __init__(self, parent = None, **kw):
	# Define the megawidget options.
	
	optiondefs = (
	    ('borderx',     20,    INITOPT),
	    ('bordery',     20,    INITOPT),
	)
	self.defineoptions(kw, optiondefs)

	# Initialise the base class (after defining the options).
	Dialog.__init__(self, parent)

	# Create the components.
	interior = self.interior()
	aliases = (
	    ('entry', 'entryfield_entry'),
	    ('label', 'entryfield_label'),
	)
	self._promptDialogEntry = self.createcomponent('entryfield',
		aliases, None,
		EntryField, (interior,))
	self._promptDialogEntry.pack(fill='x', expand=1,
		padx = self['borderx'], pady = self['bordery'])
	
        if not kw.has_key('activatecommand'):
            # Whenever this dialog is activated, set the focus to the
            # EntryField's entry widget.
            tkentry = self.component('entry')
            self.configure(activatecommand = tkentry.focus_set)

	# Check keywords and initialise options.
	self.initialiseoptions()

    # Supply aliases to some of the entry component methods.
    def insertentry(self, index, text):
	self._promptDialogEntry.insert(index, text)

    def deleteentry(self, first, last=None):
	self._promptDialogEntry.delete(first, last)

    def indexentry(self, index):
	return self._promptDialogEntry.index(index)

forwardmethods(PromptDialog, EntryField, '_promptDialogEntry')

######################################################################
### File: PmwRadioSelect.py
import types
import Tkinter


class RadioSelect(MegaWidget):
    # A collection of several buttons.  In single mode, only one
    # button may be selected.  In multiple mode, any number of buttons
    # may be selected.

    def __init__(self, parent = None, **kw):

	# Define the megawidget options.
	
	optiondefs = (
	    ('buttontype',    'button',      INITOPT),
	    ('command',       None,          None),
	    ('labelmargin',   0,             INITOPT),
	    ('labelpos',      None,          INITOPT),
	    ('orient',       'horizontal',   INITOPT),
	    ('padx',          5,             INITOPT),
	    ('pady',          5,             INITOPT),
	    ('selectmode',    'single',      INITOPT),
	)
	self.defineoptions(kw, optiondefs, dynamicGroups = ('Button',))

	# Initialise the base class (after defining the options).
	MegaWidget.__init__(self, parent)

	# Create the components.
	interior = self.interior()
	if self['labelpos'] is None:
	    self._radioSelectFrame = self._hull
	else:
	    self._radioSelectFrame = self.createcomponent('frame',
		    (), None,
		    Tkinter.Frame, (interior,))
	    self._radioSelectFrame.grid(column=2, row=2, sticky='nsew')
	    interior.grid_columnconfigure(2, weight=1)
	    interior.grid_rowconfigure(2, weight=1)

	    self.createlabel(interior)

	# Initialise instance variables.
	self._buttonList = []
	if self['selectmode'] == 'single':
	    self._singleSelect = 1
	elif self['selectmode'] == 'multiple':
	    self._singleSelect = 0
	else: 
	    raise ValueError('bad selectmode option "' + self['selectmode']
                             + '": should be single or multiple')

	if self['buttontype'] == 'button':
	    self.buttonClass = Tkinter.Button
	elif self['buttontype'] == 'radiobutton':
	    self._singleSelect = 1
	    self.var = Tkinter.StringVar()
	    self.buttonClass = Tkinter.Radiobutton
	elif self['buttontype'] == 'checkbutton':
	    self._singleSelect = 0
	    self.buttonClass = Tkinter.Checkbutton
	else:
	    raise ValueError('bad buttontype option "' + self['buttontype']
                             + '": should be button, radiobutton or checkbutton')

	if self._singleSelect:
	    self.selection = None
	else:
	    self.selection = []

	if self['orient'] not in ('horizontal', 'vertical'):
	    raise ValueError('bad orient option ' + repr(self['orient'])
                             + ': must be either \'horizontal\' or \'vertical\'')

	# Check keywords and initialise options.
	self.initialiseoptions()

    def getcurselection(self):
	if self._singleSelect:
            return self.selection
        else:
            return tuple(self.selection)

    def getvalue(self):
        return self.getcurselection()

    def setvalue(self, textOrList):
	if self._singleSelect:
            self.__setSingleValue(textOrList)
        else:
	    # Multiple selections
            oldselection = self.selection
            self.selection = textOrList
            for button in self._buttonList:
                if button in oldselection:
                    if button not in self.selection:
                        # button is currently selected but should not be
                        widget = self.component(button)
                        if self['buttontype'] == 'checkbutton':
                            widget.deselect()
                        else:  # Button
                            widget.configure(relief='raised')
                else:
                    if button in self.selection:
                        # button is not currently selected but should be
                        widget = self.component(button)
                        if self['buttontype'] == 'checkbutton':
                            widget.select()
                        else:  # Button
                            widget.configure(relief='sunken')

    def numbuttons(self):
        return len(self._buttonList)

    def index(self, index):
	# Return the integer index of the button with the given index.

	listLength = len(self._buttonList)
        if isinstance(index, int):
	    if index < listLength:
		return index
	    else:
		raise ValueError('index "%s" is out of range' % index)
	elif index is END:
	    if listLength > 0:
		return listLength - 1
	    else:
		raise ValueError('RadioSelect has no buttons')
	else:
	    for count in range(listLength):
		name = self._buttonList[count]
		if index == name:
		    return count
	    validValues = 'a name, a number or END'
	    raise ValueError('bad index "%s": must be %s' % (index, validValues))

    def button(self, buttonIndex):
	name = self._buttonList[self.index(buttonIndex)]
        return self.component(name)

    def add(self, componentName, **kw):
	if componentName in self._buttonList:
	    raise ValueError('button "%s" already exists' % componentName)

	kw['command'] = \
                lambda self=self, name=componentName: self.invoke(name)
	if not kw.has_key('text'):
	    kw['text'] = componentName

	if self['buttontype'] == 'radiobutton':
	    if not kw.has_key('anchor'):
		kw['anchor'] = 'w'
	    if not kw.has_key('variable'):
		kw['variable'] = self.var
	    if not kw.has_key('value'):
		kw['value'] = kw['text']
	elif self['buttontype'] == 'checkbutton':
	    if not kw.has_key('anchor'):
		kw['anchor'] = 'w'

	button = apply(self.createcomponent, (componentName,
		(), 'Button',
		self.buttonClass, (self._radioSelectFrame,)), kw)

	if self['orient'] == 'horizontal':
	    self._radioSelectFrame.grid_rowconfigure(0, weight=1)
	    col = len(self._buttonList)
	    button.grid(column=col, row=0, padx = self['padx'],
		    pady = self['pady'], sticky='nsew')
	    self._radioSelectFrame.grid_columnconfigure(col, weight=1)
	else:
	    self._radioSelectFrame.grid_columnconfigure(0, weight=1)
	    row = len(self._buttonList)
	    button.grid(column=0, row=row, padx = self['padx'],
		    pady = self['pady'], sticky='ew')
	    self._radioSelectFrame.grid_rowconfigure(row, weight=1)

	self._buttonList.append(componentName)
	return button

    def deleteall(self):
	for name in self._buttonList:
	    self.destroycomponent(name)
	self._buttonList = []
	if self._singleSelect:
	    self.selection = None
	else: 
	    self.selection = []

    def __setSingleValue(self, value):
            self.selection = value
            if self['buttontype'] == 'radiobutton':
                widget = self.component(value)
                widget.select()
            else:  # Button
                for button in self._buttonList:
                    widget = self.component(button)
                    if button == value:
                        widget.configure(relief='sunken')
                    else:
                        widget.configure(relief='raised')

    def invoke(self, index):
	index = self.index(index)
	name = self._buttonList[index]

	if self._singleSelect:
            self.__setSingleValue(name)
	    command = self['command']
	    if callable(command):
		return command(name)
        else:
	    # Multiple selections
	    widget = self.component(name)
	    if name in self.selection:
		if self['buttontype'] == 'checkbutton':
		    widget.deselect()
		else:
		    widget.configure(relief='raised')
		self.selection.remove(name)
		state = 0
	    else:
		if self['buttontype'] == 'checkbutton':
		    widget.select()
		else:
		    widget.configure(relief='sunken')
		self.selection.append(name)
		state = 1

	    command = self['command']
	    if callable(command):
	      return command(name, state)

######################################################################
### File: PmwScrolledCanvas.py
import Tkinter


class ScrolledCanvas(MegaWidget):
    def __init__(self, parent = None, **kw):

	# Define the megawidget options.
	
	optiondefs = (
	    ('borderframe',    0,            INITOPT),
	    ('canvasmargin',   0,            INITOPT),
	    ('hscrollmode',    'dynamic',    self._hscrollMode),
	    ('labelmargin',    0,            INITOPT),
	    ('labelpos',       None,         INITOPT),
	    ('scrollmargin',   2,            INITOPT),
	    ('usehullsize',    0,            INITOPT),
	    ('vscrollmode',    'dynamic',    self._vscrollMode),
	)
	self.defineoptions(kw, optiondefs)

	# Initialise the base class (after defining the options).
	MegaWidget.__init__(self, parent)

	# Create the components.
	self.origInterior = MegaWidget.interior(self)

	if self['usehullsize']:
	    self.origInterior.grid_propagate(0)

	if self['borderframe']:
	    # Create a frame widget to act as the border of the canvas. 
	    self._borderframe = self.createcomponent('borderframe',
		    (), None,
		    Tkinter.Frame, (self.origInterior,),
		    relief = 'sunken',
		    borderwidth = 2,
	    )
	    self._borderframe.grid(row = 2, column = 2, sticky = 'news')

	    # Create the canvas widget.
	    self._canvas = self.createcomponent('canvas',
		    (), None,
		    Tkinter.Canvas, (self._borderframe,),
		    highlightthickness = 0,
		    borderwidth = 0,
	    )
	    self._canvas.pack(fill = 'both', expand = 1)
	else:
	    # Create the canvas widget.
	    self._canvas = self.createcomponent('canvas',
		    (), None,
		    Tkinter.Canvas, (self.origInterior,),
		    relief = 'sunken',
		    borderwidth = 2,
	    )
	    self._canvas.grid(row = 2, column = 2, sticky = 'news')

	self.origInterior.grid_rowconfigure(2, weight = 1, minsize = 0)
	self.origInterior.grid_columnconfigure(2, weight = 1, minsize = 0)
	
	# Create the horizontal scrollbar
	self._horizScrollbar = self.createcomponent('horizscrollbar',
		(), 'Scrollbar',
		Tkinter.Scrollbar, (self.origInterior,),
	        orient='horizontal',
		command=self._canvas.xview
	)

	# Create the vertical scrollbar
	self._vertScrollbar = self.createcomponent('vertscrollbar',
		(), 'Scrollbar',
		Tkinter.Scrollbar, (self.origInterior,),
		orient='vertical',
		command=self._canvas.yview
	)

	self.createlabel(self.origInterior, childCols = 3, childRows = 3)

	# Initialise instance variables.
	self._horizScrollbarOn = 0
	self._vertScrollbarOn = 0
	self.scrollTimer = None
        self._scrollRecurse = 0
	self._horizScrollbarNeeded = 0
	self._vertScrollbarNeeded = 0
	self.setregionTimer = None

	# Check keywords and initialise options.
	self.initialiseoptions()

    def destroy(self):
	if self.scrollTimer is not None:
	    self.after_cancel(self.scrollTimer)
	    self.scrollTimer = None
	if self.setregionTimer is not None:
	    self.after_cancel(self.setregionTimer)
	    self.setregionTimer = None
	MegaWidget.destroy(self)

    # ======================================================================

    # Public methods.

    def interior(self):
	return self._canvas

    def resizescrollregion(self):
	if self.setregionTimer is None:
	    self.setregionTimer = self.after_idle(self._setRegion)

    # ======================================================================

    # Configuration methods.

    def _hscrollMode(self):
	# The horizontal scroll mode has been configured.

	mode = self['hscrollmode']

	if mode == 'static':
	    if not self._horizScrollbarOn:
		self._toggleHorizScrollbar()
	elif mode == 'dynamic':
	    if self._horizScrollbarNeeded != self._horizScrollbarOn:
		self._toggleHorizScrollbar()
	elif mode == 'none':
	    if self._horizScrollbarOn:
		self._toggleHorizScrollbar()
	else:
	    message = 'bad hscrollmode option "%s": should be static, dynamic, or none' % mode
	    raise ValueError(message)

        self._configureScrollCommands()

    def _vscrollMode(self):
	# The vertical scroll mode has been configured.

	mode = self['vscrollmode']

	if mode == 'static':
	    if not self._vertScrollbarOn:
		self._toggleVertScrollbar()
	elif mode == 'dynamic':
	    if self._vertScrollbarNeeded != self._vertScrollbarOn:
		self._toggleVertScrollbar()
	elif mode == 'none':
	    if self._vertScrollbarOn:
		self._toggleVertScrollbar()
	else:
	    message = 'bad vscrollmode option "%s": should be static, dynamic, or none' % mode
	    raise ValueError(message)

        self._configureScrollCommands()

    # ======================================================================

    # Private methods.

    def _configureScrollCommands(self):
        # If both scrollmodes are not dynamic we can save a lot of
        # time by not having to create an idle job to handle the
        # scroll commands.

        # Clean up previous scroll commands to prevent memory leak.
        tclCommandName = str(self._canvas.cget('xscrollcommand'))
        if tclCommandName != '':   
            self._canvas.deletecommand(tclCommandName)
        tclCommandName = str(self._canvas.cget('yscrollcommand'))
        if tclCommandName != '':   
            self._canvas.deletecommand(tclCommandName)

	if self['hscrollmode'] == self['vscrollmode'] == 'dynamic':
            self._canvas.configure(
                    xscrollcommand=self._scrollBothLater,
                    yscrollcommand=self._scrollBothLater
            )
        else:
            self._canvas.configure(
                    xscrollcommand=self._scrollXNow,
                    yscrollcommand=self._scrollYNow
            )

    def _scrollXNow(self, first, last):
        self._horizScrollbar.set(first, last)
        self._horizScrollbarNeeded = ((first, last) != ('0', '1'))

	if self['hscrollmode'] == 'dynamic':
	    if self._horizScrollbarNeeded != self._horizScrollbarOn:
		self._toggleHorizScrollbar()

    def _scrollYNow(self, first, last):
        self._vertScrollbar.set(first, last)
        self._vertScrollbarNeeded = ((first, last) != ('0', '1'))

        if self['vscrollmode'] == 'dynamic':
            if self._vertScrollbarNeeded != self._vertScrollbarOn:
                self._toggleVertScrollbar()

    def _scrollBothLater(self, first, last):
	# Called by the canvas to set the horizontal or vertical
	# scrollbar when it has scrolled or changed scrollregion.

	if self.scrollTimer is None:
	    self.scrollTimer = self.after_idle(self._scrollBothNow)

    def _scrollBothNow(self):
        # This performs the function of _scrollXNow and _scrollYNow.
        # If one is changed, the other should be updated to match.
	self.scrollTimer = None

        # Call update_idletasks to make sure that the containing frame
        # has been resized before we attempt to set the scrollbars. 
        # Otherwise the scrollbars may be mapped/unmapped continuously.
        self._scrollRecurse = self._scrollRecurse + 1
        self.update_idletasks()
        self._scrollRecurse = self._scrollRecurse - 1
        if self._scrollRecurse != 0:
            return

	xview = self._canvas.xview()
	yview = self._canvas.yview()
	self._horizScrollbar.set(xview[0], xview[1])
	self._vertScrollbar.set(yview[0], yview[1])

	self._horizScrollbarNeeded = (xview != (0.0, 1.0))
	self._vertScrollbarNeeded = (yview != (0.0, 1.0))

	# If both horizontal and vertical scrollmodes are dynamic and
	# currently only one scrollbar is mapped and both should be
	# toggled, then unmap the mapped scrollbar.  This prevents a
	# continuous mapping and unmapping of the scrollbars. 
	if (self['hscrollmode'] == self['vscrollmode'] == 'dynamic' and
		self._horizScrollbarNeeded != self._horizScrollbarOn and
		self._vertScrollbarNeeded != self._vertScrollbarOn and
		self._vertScrollbarOn != self._horizScrollbarOn):
	    if self._horizScrollbarOn:
		self._toggleHorizScrollbar()
	    else:
		self._toggleVertScrollbar()
	    return

	if self['hscrollmode'] == 'dynamic':
	    if self._horizScrollbarNeeded != self._horizScrollbarOn:
		self._toggleHorizScrollbar()

	if self['vscrollmode'] == 'dynamic':
	    if self._vertScrollbarNeeded != self._vertScrollbarOn:
		self._toggleVertScrollbar()

    def _toggleHorizScrollbar(self):

	self._horizScrollbarOn = not self._horizScrollbarOn

	interior = self.origInterior
	if self._horizScrollbarOn:
	    self._horizScrollbar.grid(row = 4, column = 2, sticky = 'news')
	    interior.grid_rowconfigure(3, minsize = self['scrollmargin'])
	else:
	    self._horizScrollbar.grid_forget()
	    interior.grid_rowconfigure(3, minsize = 0)

    def _toggleVertScrollbar(self):

	self._vertScrollbarOn = not self._vertScrollbarOn

	interior = self.origInterior
	if self._vertScrollbarOn:
	    self._vertScrollbar.grid(row = 2, column = 4, sticky = 'news')
	    interior.grid_columnconfigure(3, minsize = self['scrollmargin'])
	else:
	    self._vertScrollbar.grid_forget()
	    interior.grid_columnconfigure(3, minsize = 0)

    def _setRegion(self):
	self.setregionTimer = None

	region = self._canvas.bbox('all')
        if region is not None:
	    canvasmargin = self['canvasmargin']
	    region = (region[0] - canvasmargin, region[1] - canvasmargin,
		region[2] + canvasmargin, region[3] + canvasmargin)
	    self._canvas.configure(scrollregion = region)

    # Need to explicitly forward this to override the stupid
    # (grid_)bbox method inherited from Tkinter.Frame.Grid.
    def bbox(self, *args):
	return apply(self._canvas.bbox, args)

forwardmethods(ScrolledCanvas, Tkinter.Canvas, '_canvas')

######################################################################
### File: PmwScrolledField.py
import Tkinter


class ScrolledField(MegaWidget):
    def __init__(self, parent = None, **kw):

	# Define the megawidget options.
	
	optiondefs = (
	    ('labelmargin',   0,      INITOPT),
	    ('labelpos',      None,   INITOPT),
	    ('sticky',        'ew',   INITOPT),
	    ('text',          '',     self._text),
	)
	self.defineoptions(kw, optiondefs)

	# Initialise the base class (after defining the options).
	MegaWidget.__init__(self, parent)

	# Create the components.
	interior = self.interior()
	self._scrolledFieldEntry = self.createcomponent('entry',
		(), None,
		Tkinter.Entry, (interior,))

        # Can't always use 'disabled', since this greys out text in Tk 8.4.2
        try:
            self._scrolledFieldEntry.configure(state = 'readonly')
        except Tkinter.TclError:
            self._scrolledFieldEntry.configure(state = 'disabled')

	self._scrolledFieldEntry.grid(column=2, row=2, sticky=self['sticky'])
	interior.grid_columnconfigure(2, weight=1)
	interior.grid_rowconfigure(2, weight=1)

	self.createlabel(interior)

	# Check keywords and initialise options.
	self.initialiseoptions()

    def _text(self):
        text = self['text']
        self._scrolledFieldEntry.configure(state = 'normal')
        self._scrolledFieldEntry.delete(0, 'end')
        self._scrolledFieldEntry.insert('end', text)

        # Can't always use 'disabled', since this greys out text in Tk 8.4.2
        try:
            self._scrolledFieldEntry.configure(state = 'readonly')
        except Tkinter.TclError:
            self._scrolledFieldEntry.configure(state = 'disabled')

forwardmethods(ScrolledField, Tkinter.Entry, '_scrolledFieldEntry')

######################################################################
### File: PmwScrolledFrame.py
import string
import types
import Tkinter


class ScrolledFrame(MegaWidget):
    def __init__(self, parent = None, **kw):

	# Define the megawidget options.
	
	optiondefs = (
	    ('borderframe',    1,            INITOPT),
	    ('horizflex',      'fixed',      self._horizflex),
	    ('horizfraction',  0.05,         INITOPT),
	    ('hscrollmode',    'dynamic',    self._hscrollMode),
	    ('labelmargin',    0,            INITOPT),
	    ('labelpos',       None,         INITOPT),
	    ('scrollmargin',   2,            INITOPT),
	    ('usehullsize',    0,            INITOPT),
	    ('vertflex',       'fixed',      self._vertflex),
	    ('vertfraction',   0.05,         INITOPT),
	    ('vscrollmode',    'dynamic',    self._vscrollMode),
	)
	self.defineoptions(kw, optiondefs)

	# Initialise the base class (after defining the options).
	MegaWidget.__init__(self, parent)

	# Create the components.
	self.origInterior = MegaWidget.interior(self)

	if self['usehullsize']:
	    self.origInterior.grid_propagate(0)

	if self['borderframe']:
	    # Create a frame widget to act as the border of the clipper. 
	    self._borderframe = self.createcomponent('borderframe',
		    (), None,
		    Tkinter.Frame, (self.origInterior,),
		    relief = 'sunken',
		    borderwidth = 2,
	    )
	    self._borderframe.grid(row = 2, column = 2, sticky = 'news')

	    # Create the clipping window.
	    self._clipper = self.createcomponent('clipper',
		    (), None,
		    Tkinter.Frame, (self._borderframe,),
		    width = 400,
		    height = 300,
		    highlightthickness = 0,
		    borderwidth = 0,
	    )
	    self._clipper.pack(fill = 'both', expand = 1)
	else:
	    # Create the clipping window.
	    self._clipper = self.createcomponent('clipper',
		    (), None,
		    Tkinter.Frame, (self.origInterior,),
		    width = 400,
		    height = 300,
		    relief = 'sunken',
		    borderwidth = 2,
	    )
	    self._clipper.grid(row = 2, column = 2, sticky = 'news')

	self.origInterior.grid_rowconfigure(2, weight = 1, minsize = 0)
	self.origInterior.grid_columnconfigure(2, weight = 1, minsize = 0)
	
	# Create the horizontal scrollbar
	self._horizScrollbar = self.createcomponent('horizscrollbar',
		(), 'Scrollbar',
		Tkinter.Scrollbar, (self.origInterior,),
	        orient='horizontal',
		command=self.xview
	)

	# Create the vertical scrollbar
	self._vertScrollbar = self.createcomponent('vertscrollbar',
		(), 'Scrollbar',
		Tkinter.Scrollbar, (self.origInterior,),
		orient='vertical',
		command=self.yview
	)

	self.createlabel(self.origInterior, childCols = 3, childRows = 3)

	# Initialise instance variables.
	self._horizScrollbarOn = 0
	self._vertScrollbarOn = 0
	self.scrollTimer = None
	self._scrollRecurse = 0
	self._horizScrollbarNeeded = 0
	self._vertScrollbarNeeded = 0
	self.startX = 0
	self.startY = 0
	self._flexoptions = ('fixed', 'expand', 'shrink', 'elastic')

	# Create a frame in the clipper to contain the widgets to be
	# scrolled.
	self._frame = self.createcomponent('frame',
		(), None,
		Tkinter.Frame, (self._clipper,)
	)

	# Whenever the clipping window or scrolled frame change size,
	# update the scrollbars.
	self._frame.bind('<Configure>', self._reposition)
	self._clipper.bind('<Configure>', self._reposition)

        # Work around a bug in Tk where the value returned by the
        # scrollbar get() method is (0.0, 0.0, 0.0, 0.0) rather than
        # the expected 2-tuple.  This occurs if xview() is called soon
        # after the ScrolledFrame has been created.
        self._horizScrollbar.set(0.0, 1.0)
        self._vertScrollbar.set(0.0, 1.0)

	# Check keywords and initialise options.
	self.initialiseoptions()

    def destroy(self):
	if self.scrollTimer is not None:
	    self.after_cancel(self.scrollTimer)
	    self.scrollTimer = None
	MegaWidget.destroy(self)

    # ======================================================================

    # Public methods.

    def interior(self):
	return self._frame

    # Set timer to call real reposition method, so that it is not
    # called multiple times when many things are reconfigured at the
    # same time.
    def reposition(self):
	if self.scrollTimer is None:
	    self.scrollTimer = self.after_idle(self._scrollBothNow)

    # Called when the user clicks in the horizontal scrollbar. 
    # Calculates new position of frame then calls reposition() to
    # update the frame and the scrollbar.
    def xview(self, mode = None, value = None, units = None):

        if isinstance(value, basestring):
            value = string.atof(value)
        if mode is None:
            return self._horizScrollbar.get()
	elif mode == 'moveto':
	    frameWidth = self._frame.winfo_reqwidth()
	    self.startX = value * float(frameWidth)
	else: # mode == 'scroll'
	    clipperWidth = self._clipper.winfo_width()
	    if units == 'units':
		jump = int(clipperWidth * self['horizfraction'])
	    else:
		jump = clipperWidth
            self.startX = self.startX + value * jump

	self.reposition()

    # Called when the user clicks in the vertical scrollbar. 
    # Calculates new position of frame then calls reposition() to
    # update the frame and the scrollbar.
    def yview(self, mode = None, value = None, units = None):

        if isinstance(value, basestring):
            value = string.atof(value)
        if mode is None:
            return self._vertScrollbar.get()
	elif mode == 'moveto':
	    frameHeight = self._frame.winfo_reqheight()
	    self.startY = value * float(frameHeight)
	else: # mode == 'scroll'
	    clipperHeight = self._clipper.winfo_height()
	    if units == 'units':
		jump = int(clipperHeight * self['vertfraction'])
	    else:
		jump = clipperHeight
            self.startY = self.startY + value * jump

	self.reposition()

    # ======================================================================

    # Configuration methods.

    def _hscrollMode(self):
	# The horizontal scroll mode has been configured.

	mode = self['hscrollmode']

	if mode == 'static':
	    if not self._horizScrollbarOn:
		self._toggleHorizScrollbar()
	elif mode == 'dynamic':
	    if self._horizScrollbarNeeded != self._horizScrollbarOn:
		self._toggleHorizScrollbar()
	elif mode == 'none':
	    if self._horizScrollbarOn:
		self._toggleHorizScrollbar()
	else:
	    message = 'bad hscrollmode option "%s": should be static, dynamic, or none' % mode
	    raise ValueError(message)

    def _vscrollMode(self):
	# The vertical scroll mode has been configured.

	mode = self['vscrollmode']

	if mode == 'static':
	    if not self._vertScrollbarOn:
		self._toggleVertScrollbar()
	elif mode == 'dynamic':
	    if self._vertScrollbarNeeded != self._vertScrollbarOn:
		self._toggleVertScrollbar()
	elif mode == 'none':
	    if self._vertScrollbarOn:
		self._toggleVertScrollbar()
	else:
	    message = 'bad vscrollmode option "%s": should be static, dynamic, or none' % mode
	    raise ValueError(message)

    def _horizflex(self):
	# The horizontal flex mode has been configured.

	flex = self['horizflex']

	if flex not in self._flexoptions:
	    message = 'bad horizflex option "%s": should be one of %s' % \
		    (flex, str(self._flexoptions))
	    raise ValueError(message)

	self.reposition()

    def _vertflex(self):
	# The vertical flex mode has been configured.

	flex = self['vertflex']

	if flex not in self._flexoptions:
	    message = 'bad vertflex option "%s": should be one of %s' % \
		    (flex, str(self._flexoptions))
	    raise ValueError(message)

	self.reposition()

    # ======================================================================

    # Private methods.

    def _reposition(self, event):
	self.reposition()

    def _getxview(self):

	# Horizontal dimension.
	clipperWidth = self._clipper.winfo_width()
	frameWidth = self._frame.winfo_reqwidth()
	if frameWidth <= clipperWidth:
	    # The scrolled frame is smaller than the clipping window.

	    self.startX = 0
	    endScrollX = 1.0

	    if self['horizflex'] in ('expand', 'elastic'):
		relwidth = 1
	    else:
		relwidth = ''
	else:
	    # The scrolled frame is larger than the clipping window.

	    if self['horizflex'] in ('shrink', 'elastic'):
		self.startX = 0
		endScrollX = 1.0
		relwidth = 1
	    else:
		if self.startX + clipperWidth > frameWidth:
		    self.startX = frameWidth - clipperWidth
		    endScrollX = 1.0
		else:
		    if self.startX < 0:
			self.startX = 0
		    endScrollX = (self.startX + clipperWidth) / float(frameWidth)
		relwidth = ''

	# Position frame relative to clipper.
	self._frame.place(x = -self.startX, relwidth = relwidth)
	return (self.startX / float(frameWidth), endScrollX)

    def _getyview(self):

	# Vertical dimension.
	clipperHeight = self._clipper.winfo_height()
	frameHeight = self._frame.winfo_reqheight()
	if frameHeight <= clipperHeight:
	    # The scrolled frame is smaller than the clipping window.

	    self.startY = 0
	    endScrollY = 1.0

	    if self['vertflex'] in ('expand', 'elastic'):
		relheight = 1
	    else:
		relheight = ''
	else:
	    # The scrolled frame is larger than the clipping window.

	    if self['vertflex'] in ('shrink', 'elastic'):
		self.startY = 0
		endScrollY = 1.0
		relheight = 1
	    else:
		if self.startY + clipperHeight > frameHeight:
		    self.startY = frameHeight - clipperHeight
		    endScrollY = 1.0
		else:
		    if self.startY < 0:
			self.startY = 0
		    endScrollY = (self.startY + clipperHeight) / float(frameHeight)
		relheight = ''

	# Position frame relative to clipper.
	self._frame.place(y = -self.startY, relheight = relheight)
	return (self.startY / float(frameHeight), endScrollY)

    # According to the relative geometries of the frame and the
    # clipper, reposition the frame within the clipper and reset the
    # scrollbars.
    def _scrollBothNow(self):
	self.scrollTimer = None

        # Call update_idletasks to make sure that the containing frame
        # has been resized before we attempt to set the scrollbars. 
        # Otherwise the scrollbars may be mapped/unmapped continuously.
        self._scrollRecurse = self._scrollRecurse + 1
        self.update_idletasks()
        self._scrollRecurse = self._scrollRecurse - 1
        if self._scrollRecurse != 0:
            return

	xview = self._getxview()
	yview = self._getyview()
	self._horizScrollbar.set(xview[0], xview[1])
	self._vertScrollbar.set(yview[0], yview[1])

	self._horizScrollbarNeeded = (xview != (0.0, 1.0))
	self._vertScrollbarNeeded = (yview != (0.0, 1.0))

	# If both horizontal and vertical scrollmodes are dynamic and
	# currently only one scrollbar is mapped and both should be
	# toggled, then unmap the mapped scrollbar.  This prevents a
	# continuous mapping and unmapping of the scrollbars. 
	if (self['hscrollmode'] == self['vscrollmode'] == 'dynamic' and
		self._horizScrollbarNeeded != self._horizScrollbarOn and
		self._vertScrollbarNeeded != self._vertScrollbarOn and
		self._vertScrollbarOn != self._horizScrollbarOn):
	    if self._horizScrollbarOn:
		self._toggleHorizScrollbar()
	    else:
		self._toggleVertScrollbar()
	    return

	if self['hscrollmode'] == 'dynamic':
	    if self._horizScrollbarNeeded != self._horizScrollbarOn:
		self._toggleHorizScrollbar()

	if self['vscrollmode'] == 'dynamic':
	    if self._vertScrollbarNeeded != self._vertScrollbarOn:
		self._toggleVertScrollbar()

    def _toggleHorizScrollbar(self):

	self._horizScrollbarOn = not self._horizScrollbarOn

	interior = self.origInterior
	if self._horizScrollbarOn:
	    self._horizScrollbar.grid(row = 4, column = 2, sticky = 'news')
	    interior.grid_rowconfigure(3, minsize = self['scrollmargin'])
	else:
	    self._horizScrollbar.grid_forget()
	    interior.grid_rowconfigure(3, minsize = 0)

    def _toggleVertScrollbar(self):

	self._vertScrollbarOn = not self._vertScrollbarOn

	interior = self.origInterior
	if self._vertScrollbarOn:
	    self._vertScrollbar.grid(row = 2, column = 4, sticky = 'news')
	    interior.grid_columnconfigure(3, minsize = self['scrollmargin'])
	else:
	    self._vertScrollbar.grid_forget()
	    interior.grid_columnconfigure(3, minsize = 0)

######################################################################
### File: PmwScrolledListBox.py
# Based on iwidgets2.2.0/scrolledlistbox.itk code.

import types
import Tkinter


class ScrolledListBox(MegaWidget):
    _classBindingsDefinedFor = 0

    def __init__(self, parent = None, **kw):

	# Define the megawidget options.
	
	optiondefs = (
	    ('dblclickcommand',    None,            None),
	    ('hscrollmode',        'dynamic',       self._hscrollMode),
	    ('items',              (),              INITOPT),
	    ('labelmargin',        0,               INITOPT),
	    ('labelpos',           None,            INITOPT),
	    ('scrollmargin',       2,               INITOPT),
	    ('selectioncommand',   None,            None),
	    ('usehullsize',        0,               INITOPT),
	    ('vscrollmode',        'dynamic',       self._vscrollMode),
	)
	self.defineoptions(kw, optiondefs)

	# Initialise the base class (after defining the options).
	MegaWidget.__init__(self, parent)

	# Create the components.
	interior = self.interior()

	if self['usehullsize']:
	    interior.grid_propagate(0)

	# Create the listbox widget.
	self._listbox = self.createcomponent('listbox',
		(), None,
		Tkinter.Listbox, (interior,))
	self._listbox.grid(row = 2, column = 2, sticky = 'news')
	interior.grid_rowconfigure(2, weight = 1, minsize = 0)
	interior.grid_columnconfigure(2, weight = 1, minsize = 0)

	# Create the horizontal scrollbar
	self._horizScrollbar = self.createcomponent('horizscrollbar',
		(), 'Scrollbar',
		Tkinter.Scrollbar, (interior,),
	        orient='horizontal',
		command=self._listbox.xview
	)

	# Create the vertical scrollbar
	self._vertScrollbar = self.createcomponent('vertscrollbar',
		(), 'Scrollbar',
		Tkinter.Scrollbar, (interior,),
		orient='vertical',
		command=self._listbox.yview
	)

	self.createlabel(interior, childCols = 3, childRows = 3)

	# Add the items specified by the initialisation option.
	items = self['items']
        if not isinstance(items, tuple):
	    items = tuple(items)
	if len(items) > 0:
	    apply(self._listbox.insert, ('end',) + items)

	_registerScrolledList(self._listbox, self)

        # Establish the special class bindings if not already done.
        # Also create bindings if the Tkinter default interpreter has
        # changed.  Use Tkinter._default_root to create class
        # bindings, so that a reference to root is created by
        # bind_class rather than a reference to self, which would
        # prevent object cleanup.
        theTag = 'ScrolledListBoxTag'
        if ScrolledListBox._classBindingsDefinedFor != Tkinter._default_root:
            root  = Tkinter._default_root
	    	    
            def doubleEvent(event):
                _handleEvent(event, 'double')
            def keyEvent(event):
                _handleEvent(event, 'key')
            def releaseEvent(event):
                _handleEvent(event, 'release')

            # Bind space and return keys and button 1 to the selectioncommand.
            root.bind_class(theTag, '<Key-space>', keyEvent)
            root.bind_class(theTag, '<Key-Return>', keyEvent)
            root.bind_class(theTag, '<ButtonRelease-1>', releaseEvent)

            # Bind double button 1 click to the dblclickcommand.
            root.bind_class(theTag, '<Double-ButtonRelease-1>', doubleEvent)

	    ScrolledListBox._classBindingsDefinedFor = root

	bindtags = self._listbox.bindtags()
	self._listbox.bindtags(bindtags + (theTag,))

	# Initialise instance variables.
	self._horizScrollbarOn = 0
	self._vertScrollbarOn = 0
	self.scrollTimer = None
        self._scrollRecurse = 0
	self._horizScrollbarNeeded = 0
	self._vertScrollbarNeeded = 0

	# Check keywords and initialise options.
	self.initialiseoptions()

    def destroy(self):
	if self.scrollTimer is not None:
	    self.after_cancel(self.scrollTimer)
	    self.scrollTimer = None
	_deregisterScrolledList(self._listbox)
	MegaWidget.destroy(self)

    # ======================================================================

    # Public methods.

    def clear(self):
	self.setlist(())

    def getcurselection(self):
	rtn = []
	for sel in self.curselection():
	    rtn.append(self._listbox.get(sel))
	return tuple(rtn)

    def getvalue(self):
        return self.getcurselection()

    def setvalue(self, textOrList):
        self._listbox.selection_clear(0, 'end')
        listitems = list(self._listbox.get(0, 'end'))
        if isinstance(textOrList, basestring):
            if textOrList in listitems:
                self._listbox.selection_set(listitems.index(textOrList))
            else:
                raise ValueError('no such item "%s"' % textOrList)
        else:
            for item in textOrList:
                if item in listitems:
                    self._listbox.selection_set(listitems.index(item))
                else:
                    raise ValueError('no such item "%s"' % item)

    def setlist(self, items):
        self._listbox.delete(0, 'end')
	if len(items) > 0:
            if not isinstance(items, tuple):
		items = tuple(items)
	    apply(self._listbox.insert, (0,) + items)

    # Override Tkinter.Listbox get method, so that if it is called with
    # no arguments, return all list elements (consistent with other widgets).
    def get(self, first=None, last=None):
	if first is None:
	    return self._listbox.get(0, 'end')
	else:
	    return self._listbox.get(first, last)

    # ======================================================================

    # Configuration methods.

    def _hscrollMode(self):
	# The horizontal scroll mode has been configured.

	mode = self['hscrollmode']

	if mode == 'static':
	    if not self._horizScrollbarOn:
		self._toggleHorizScrollbar()
	elif mode == 'dynamic':
	    if self._horizScrollbarNeeded != self._horizScrollbarOn:
		self._toggleHorizScrollbar()
	elif mode == 'none':
	    if self._horizScrollbarOn:
		self._toggleHorizScrollbar()
	else:
	    message = 'bad hscrollmode option "%s": should be static, dynamic, or none' % mode
	    raise ValueError(message)

        self._configureScrollCommands()

    def _vscrollMode(self):
	# The vertical scroll mode has been configured.

	mode = self['vscrollmode']

	if mode == 'static':
	    if not self._vertScrollbarOn:
		self._toggleVertScrollbar()
	elif mode == 'dynamic':
	    if self._vertScrollbarNeeded != self._vertScrollbarOn:
		self._toggleVertScrollbar()
	elif mode == 'none':
	    if self._vertScrollbarOn:
		self._toggleVertScrollbar()
	else:
	    message = 'bad vscrollmode option "%s": should be static, dynamic, or none' % mode
	    raise ValueError(message)

        self._configureScrollCommands()

    # ======================================================================

    # Private methods.

    def _configureScrollCommands(self):
        # If both scrollmodes are not dynamic we can save a lot of
        # time by not having to create an idle job to handle the
        # scroll commands.

        # Clean up previous scroll commands to prevent memory leak.
        tclCommandName = str(self._listbox.cget('xscrollcommand'))
        if tclCommandName != '':   
            self._listbox.deletecommand(tclCommandName)
        tclCommandName = str(self._listbox.cget('yscrollcommand'))
        if tclCommandName != '':   
            self._listbox.deletecommand(tclCommandName)

	if self['hscrollmode'] == self['vscrollmode'] == 'dynamic':
            self._listbox.configure(
                    xscrollcommand=self._scrollBothLater,
                    yscrollcommand=self._scrollBothLater
            )
        else:
            self._listbox.configure(
                    xscrollcommand=self._scrollXNow,
                    yscrollcommand=self._scrollYNow
            )

    def _scrollXNow(self, first, last):
        self._horizScrollbar.set(first, last)
        self._horizScrollbarNeeded = ((first, last) != ('0', '1'))

	if self['hscrollmode'] == 'dynamic':
	    if self._horizScrollbarNeeded != self._horizScrollbarOn:
		self._toggleHorizScrollbar()

    def _scrollYNow(self, first, last):
        self._vertScrollbar.set(first, last)
        self._vertScrollbarNeeded = ((first, last) != ('0', '1'))

        if self['vscrollmode'] == 'dynamic':
            if self._vertScrollbarNeeded != self._vertScrollbarOn:
                self._toggleVertScrollbar()

    def _scrollBothLater(self, first, last):
	# Called by the listbox to set the horizontal or vertical
	# scrollbar when it has scrolled or changed size or contents.

	if self.scrollTimer is None:
	    self.scrollTimer = self.after_idle(self._scrollBothNow)

    def _scrollBothNow(self):
        # This performs the function of _scrollXNow and _scrollYNow.
        # If one is changed, the other should be updated to match.
	self.scrollTimer = None

        # Call update_idletasks to make sure that the containing frame
        # has been resized before we attempt to set the scrollbars. 
        # Otherwise the scrollbars may be mapped/unmapped continuously.
        self._scrollRecurse = self._scrollRecurse + 1
        self.update_idletasks()
        self._scrollRecurse = self._scrollRecurse - 1
        if self._scrollRecurse != 0:
            return

	xview = self._listbox.xview()
	yview = self._listbox.yview()
	self._horizScrollbar.set(xview[0], xview[1])
	self._vertScrollbar.set(yview[0], yview[1])

	self._horizScrollbarNeeded = (xview != (0.0, 1.0))
	self._vertScrollbarNeeded = (yview != (0.0, 1.0))

	# If both horizontal and vertical scrollmodes are dynamic and
	# currently only one scrollbar is mapped and both should be
	# toggled, then unmap the mapped scrollbar.  This prevents a
	# continuous mapping and unmapping of the scrollbars. 
	if (self['hscrollmode'] == self['vscrollmode'] == 'dynamic' and
		self._horizScrollbarNeeded != self._horizScrollbarOn and
		self._vertScrollbarNeeded != self._vertScrollbarOn and
		self._vertScrollbarOn != self._horizScrollbarOn):
	    if self._horizScrollbarOn:
		self._toggleHorizScrollbar()
	    else:
		self._toggleVertScrollbar()
	    return

	if self['hscrollmode'] == 'dynamic':
	    if self._horizScrollbarNeeded != self._horizScrollbarOn:
		self._toggleHorizScrollbar()

	if self['vscrollmode'] == 'dynamic':
	    if self._vertScrollbarNeeded != self._vertScrollbarOn:
		self._toggleVertScrollbar()

    def _toggleHorizScrollbar(self):

	self._horizScrollbarOn = not self._horizScrollbarOn

	interior = self.interior()
	if self._horizScrollbarOn:
	    self._horizScrollbar.grid(row = 4, column = 2, sticky = 'news')
	    interior.grid_rowconfigure(3, minsize = self['scrollmargin'])
	else:
	    self._horizScrollbar.grid_forget()
	    interior.grid_rowconfigure(3, minsize = 0)

    def _toggleVertScrollbar(self):

	self._vertScrollbarOn = not self._vertScrollbarOn

	interior = self.interior()
	if self._vertScrollbarOn:
	    self._vertScrollbar.grid(row = 2, column = 4, sticky = 'news')
	    interior.grid_columnconfigure(3, minsize = self['scrollmargin'])
	else:
	    self._vertScrollbar.grid_forget()
	    interior.grid_columnconfigure(3, minsize = 0)

    def _handleEvent(self, event, eventType):
        if eventType == 'double':
            command = self['dblclickcommand']
        elif eventType == 'key':
            command = self['selectioncommand']
        else: #eventType == 'release'
            # Do not execute the command if the mouse was released
            # outside the listbox.
            if (event.x < 0 or self._listbox.winfo_width() <= event.x or
                    event.y < 0 or self._listbox.winfo_height() <= event.y):
                return

            command = self['selectioncommand']

        if callable(command):
            command()

    # Need to explicitly forward this to override the stupid
    # (grid_)size method inherited from Tkinter.Frame.Grid.
    def size(self):
	return self._listbox.size()

    # Need to explicitly forward this to override the stupid
    # (grid_)bbox method inherited from Tkinter.Frame.Grid.
    def bbox(self, index):
	return self._listbox.bbox(index)

forwardmethods(ScrolledListBox, Tkinter.Listbox, '_listbox')

# ======================================================================

_listboxCache = {}

def _registerScrolledList(listbox, scrolledList):
    # Register an ScrolledList widget for a Listbox widget

    _listboxCache[listbox] = scrolledList

def _deregisterScrolledList(listbox):
    # Deregister a Listbox widget
    del _listboxCache[listbox]

def _handleEvent(event, eventType):
    # Forward events for a Listbox to it's ScrolledListBox

    # A binding earlier in the bindtags list may have destroyed the
    # megawidget, so need to check.
    if _listboxCache.has_key(event.widget):
        _listboxCache[event.widget]._handleEvent(event, eventType)

######################################################################
### File: PmwScrolledText.py
# Based on iwidgets2.2.0/scrolledtext.itk code.   

import Tkinter


class ScrolledText(MegaWidget):
    def __init__(self, parent = None, **kw):

	# Define the megawidget options.
	
	optiondefs = (
	    ('borderframe',    0,            INITOPT),
	    ('columnheader',   0,            INITOPT),
	    ('hscrollmode',    'dynamic',    self._hscrollMode),
	    ('labelmargin',    0,            INITOPT),
	    ('labelpos',       None,         INITOPT),
	    ('rowcolumnheader',0,            INITOPT),
	    ('rowheader',      0,            INITOPT),
	    ('scrollmargin',   2,            INITOPT),
	    ('usehullsize',    0,            INITOPT),
	    ('vscrollmode',    'dynamic',    self._vscrollMode),
	)
	self.defineoptions(kw, optiondefs)

	# Initialise the base class (after defining the options).
	MegaWidget.__init__(self, parent)

	# Create the components.
	interior = self.interior()

	if self['usehullsize']:
	    interior.grid_propagate(0)

	if self['borderframe']:
	    # Create a frame widget to act as the border of the text 
	    # widget.  Later, pack the text widget so that it fills
	    # the frame.  This avoids a problem in Tk, where window
	    # items in a text widget may overlap the border of the
	    # text widget.
	    self._borderframe = self.createcomponent('borderframe',
		    (), None,
		    Tkinter.Frame, (interior,),
		    relief = 'sunken',
		    borderwidth = 2,
	    )
	    self._borderframe.grid(row = 4, column = 4, sticky = 'news')

	    # Create the text widget.
	    self._textbox = self.createcomponent('text',
		    (), None,
		    Tkinter.Text, (self._borderframe,),
		    highlightthickness = 0,
		    borderwidth = 0,
	    )
	    self._textbox.pack(fill = 'both', expand = 1)

            bw = self._borderframe.cget('borderwidth'),
            ht = self._borderframe.cget('highlightthickness'),
	else:
	    # Create the text widget.
	    self._textbox = self.createcomponent('text',
		    (), None,
		    Tkinter.Text, (interior,),
	    )
	    self._textbox.grid(row = 4, column = 4, sticky = 'news')

            bw = self._textbox.cget('borderwidth'),
            ht = self._textbox.cget('highlightthickness'),

        # Create the header text widgets
        if self['columnheader']:
            self._columnheader = self.createcomponent('columnheader',
                    (), 'Header',
                    Tkinter.Text, (interior,),
                    height=1,
                    wrap='none',
                    borderwidth = bw,
                    highlightthickness = ht,
            )
            self._columnheader.grid(row = 2, column = 4, sticky = 'ew')
            self._columnheader.configure(
                    xscrollcommand = self._columnheaderscrolled)

        if self['rowheader']:
            self._rowheader = self.createcomponent('rowheader',
                    (), 'Header',
                    Tkinter.Text, (interior,),
                    wrap='none',
                    borderwidth = bw,
                    highlightthickness = ht,
            )
            self._rowheader.grid(row = 4, column = 2, sticky = 'ns')
            self._rowheader.configure(
                    yscrollcommand = self._rowheaderscrolled)

        if self['rowcolumnheader']:
            self._rowcolumnheader = self.createcomponent('rowcolumnheader',
                    (), 'Header',
                    Tkinter.Text, (interior,),
                    height=1,
                    wrap='none',
                    borderwidth = bw,
                    highlightthickness = ht,
            )
            self._rowcolumnheader.grid(row = 2, column = 2, sticky = 'nsew')

	interior.grid_rowconfigure(4, weight = 1, minsize = 0)
	interior.grid_columnconfigure(4, weight = 1, minsize = 0)

	# Create the horizontal scrollbar
	self._horizScrollbar = self.createcomponent('horizscrollbar',
		(), 'Scrollbar',
		Tkinter.Scrollbar, (interior,),
	        orient='horizontal',
		command=self._textbox.xview
	)

	# Create the vertical scrollbar
	self._vertScrollbar = self.createcomponent('vertscrollbar',
		(), 'Scrollbar',
		Tkinter.Scrollbar, (interior,),
		orient='vertical',
		command=self._textbox.yview
	)

	self.createlabel(interior, childCols = 5, childRows = 5)

	# Initialise instance variables.
	self._horizScrollbarOn = 0
	self._vertScrollbarOn = 0
	self.scrollTimer = None
        self._scrollRecurse = 0
	self._horizScrollbarNeeded = 0
	self._vertScrollbarNeeded = 0
	self._textWidth = None

        # These four variables avoid an infinite loop caused by the
        # row or column header's scrollcommand causing the main text
        # widget's scrollcommand to be called and vice versa.
	self._textboxLastX = None
	self._textboxLastY = None
	self._columnheaderLastX = None
	self._rowheaderLastY = None

	# Check keywords and initialise options.
	self.initialiseoptions()

    def destroy(self):
	if self.scrollTimer is not None:
	    self.after_cancel(self.scrollTimer)
	    self.scrollTimer = None
	MegaWidget.destroy(self)

    # ======================================================================

    # Public methods.

    def clear(self):
	self.settext('')

    def importfile(self, fileName, where = 'end'):
	file = open(fileName, 'r')
	self._textbox.insert(where, file.read())
	file.close()

    def exportfile(self, fileName):
	file = open(fileName, 'w')
	file.write(self._textbox.get('1.0', 'end'))
	file.close()

    def settext(self, text):
	disabled = (str(self._textbox.cget('state')) == 'disabled')
	if disabled:
	    self._textbox.configure(state='normal')
	self._textbox.delete('0.0', 'end')
	self._textbox.insert('end', text)
	if disabled:
	    self._textbox.configure(state='disabled')

    # Override Tkinter.Text get method, so that if it is called with
    # no arguments, return all text (consistent with other widgets).
    def get(self, first=None, last=None):
	if first is None:
	    return self._textbox.get('1.0', 'end')
	else:
	    return self._textbox.get(first, last)

    def getvalue(self):
        return self.get()

    def setvalue(self, text):
        return self.settext(text)

    def appendtext(self, text):
        oldTop, oldBottom = self._textbox.yview()
     
        disabled = (str(self._textbox.cget('state')) == 'disabled')
        if disabled:
            self._textbox.configure(state='normal')
        self._textbox.insert('end', text)
        if disabled:
            self._textbox.configure(state='disabled')
     
        if oldBottom == 1.0:
            self._textbox.yview('moveto', 1.0)

    # ======================================================================

    # Configuration methods.

    def _hscrollMode(self):
	# The horizontal scroll mode has been configured.

	mode = self['hscrollmode']

	if mode == 'static':
	    if not self._horizScrollbarOn:
		self._toggleHorizScrollbar()
	elif mode == 'dynamic':
	    if self._horizScrollbarNeeded != self._horizScrollbarOn:
		self._toggleHorizScrollbar()
	elif mode == 'none':
	    if self._horizScrollbarOn:
		self._toggleHorizScrollbar()
	else:
	    message = 'bad hscrollmode option "%s": should be static, dynamic, or none' % mode
	    raise ValueError(message)

        self._configureScrollCommands()

    def _vscrollMode(self):
	# The vertical scroll mode has been configured.

	mode = self['vscrollmode']

	if mode == 'static':
	    if not self._vertScrollbarOn:
		self._toggleVertScrollbar()
	elif mode == 'dynamic':
	    if self._vertScrollbarNeeded != self._vertScrollbarOn:
		self._toggleVertScrollbar()
	elif mode == 'none':
	    if self._vertScrollbarOn:
		self._toggleVertScrollbar()
	else:
	    message = 'bad vscrollmode option "%s": should be static, dynamic, or none' % mode
	    raise ValueError(message)

        self._configureScrollCommands()

    # ======================================================================

    # Private methods.

    def _configureScrollCommands(self):
        # If both scrollmodes are not dynamic we can save a lot of
        # time by not having to create an idle job to handle the
        # scroll commands.

        # Clean up previous scroll commands to prevent memory leak.
        tclCommandName = str(self._textbox.cget('xscrollcommand'))
        if tclCommandName != '':   
            self._textbox.deletecommand(tclCommandName)
        tclCommandName = str(self._textbox.cget('yscrollcommand'))
        if tclCommandName != '':   
            self._textbox.deletecommand(tclCommandName)

	if self['hscrollmode'] == self['vscrollmode'] == 'dynamic':
            self._textbox.configure(
                    xscrollcommand=self._scrollBothLater,
                    yscrollcommand=self._scrollBothLater
            )
        else:
            self._textbox.configure(
                    xscrollcommand=self._scrollXNow,
                    yscrollcommand=self._scrollYNow
            )

    def _scrollXNow(self, first, last):
        self._horizScrollbar.set(first, last)
        self._horizScrollbarNeeded = ((first, last) != ('0', '1'))

        # This code is the same as in _scrollBothNow.  Keep it that way.
        if self['hscrollmode'] == 'dynamic':
            currentWidth = self._textbox.winfo_width()
            if self._horizScrollbarNeeded != self._horizScrollbarOn:
                if self._horizScrollbarNeeded or \
                        self._textWidth != currentWidth:
                    self._toggleHorizScrollbar()
            self._textWidth = currentWidth

        if self['columnheader']:
	    if self._columnheaderLastX != first:
		self._columnheaderLastX = first
		self._columnheader.xview('moveto', first)

    def _scrollYNow(self, first, last):
        if first == '0' and last == '0':
            return
        self._vertScrollbar.set(first, last)
        self._vertScrollbarNeeded = ((first, last) != ('0', '1'))

        if self['vscrollmode'] == 'dynamic':
            if self._vertScrollbarNeeded != self._vertScrollbarOn:
                self._toggleVertScrollbar()

        if self['rowheader']:
	    if self._rowheaderLastY != first:
		self._rowheaderLastY = first
		self._rowheader.yview('moveto', first)

    def _scrollBothLater(self, first, last):
	# Called by the text widget to set the horizontal or vertical
	# scrollbar when it has scrolled or changed size or contents.

	if self.scrollTimer is None:
	    self.scrollTimer = self.after_idle(self._scrollBothNow)

    def _scrollBothNow(self):
        # This performs the function of _scrollXNow and _scrollYNow.
        # If one is changed, the other should be updated to match.
	self.scrollTimer = None

        # Call update_idletasks to make sure that the containing frame
        # has been resized before we attempt to set the scrollbars. 
        # Otherwise the scrollbars may be mapped/unmapped continuously.
        self._scrollRecurse = self._scrollRecurse + 1
        self.update_idletasks()
        self._scrollRecurse = self._scrollRecurse - 1
        if self._scrollRecurse != 0:
            return

	xview = self._textbox.xview()
	yview = self._textbox.yview()

	# The text widget returns a yview of (0.0, 0.0) just after it
	# has been created. Ignore this.
	if yview == (0.0, 0.0):
	    return

        if self['columnheader']:
	    if self._columnheaderLastX != xview[0]:
		self._columnheaderLastX = xview[0]
		self._columnheader.xview('moveto', xview[0])
        if self['rowheader']:
	    if self._rowheaderLastY != yview[0]:
		self._rowheaderLastY = yview[0]
		self._rowheader.yview('moveto', yview[0])

	self._horizScrollbar.set(xview[0], xview[1])
	self._vertScrollbar.set(yview[0], yview[1])

	self._horizScrollbarNeeded = (xview != (0.0, 1.0))
	self._vertScrollbarNeeded = (yview != (0.0, 1.0))

	# If both horizontal and vertical scrollmodes are dynamic and
	# currently only one scrollbar is mapped and both should be
	# toggled, then unmap the mapped scrollbar.  This prevents a
	# continuous mapping and unmapping of the scrollbars. 
	if (self['hscrollmode'] == self['vscrollmode'] == 'dynamic' and
		self._horizScrollbarNeeded != self._horizScrollbarOn and
		self._vertScrollbarNeeded != self._vertScrollbarOn and
		self._vertScrollbarOn != self._horizScrollbarOn):
	    if self._horizScrollbarOn:
		self._toggleHorizScrollbar()
	    else:
		self._toggleVertScrollbar()
	    return

	if self['hscrollmode'] == 'dynamic':

	    # The following test is done to prevent continuous
	    # mapping and unmapping of the horizontal scrollbar. 
	    # This may occur when some event (scrolling, resizing
	    # or text changes) modifies the displayed text such
	    # that the bottom line in the window is the longest
	    # line displayed.  If this causes the horizontal
	    # scrollbar to be mapped, the scrollbar may "cover up"
	    # the bottom line, which would mean that the scrollbar
	    # is no longer required.  If the scrollbar is then
	    # unmapped, the bottom line will then become visible
	    # again, which would cause the scrollbar to be mapped
	    # again, and so on...
	    #
	    # The idea is that, if the width of the text widget
	    # has not changed and the scrollbar is currently
	    # mapped, then do not unmap the scrollbar even if it
	    # is no longer required.  This means that, during
	    # normal scrolling of the text, once the horizontal
	    # scrollbar has been mapped it will not be unmapped
	    # (until the width of the text widget changes).

	    currentWidth = self._textbox.winfo_width()
	    if self._horizScrollbarNeeded != self._horizScrollbarOn:
		if self._horizScrollbarNeeded or \
			self._textWidth != currentWidth:
		    self._toggleHorizScrollbar()
	    self._textWidth = currentWidth

	if self['vscrollmode'] == 'dynamic':
	    if self._vertScrollbarNeeded != self._vertScrollbarOn:
		self._toggleVertScrollbar()

    def _columnheaderscrolled(self, first, last):
	if self._textboxLastX != first:
	    self._textboxLastX = first
	    self._textbox.xview('moveto', first)

    def _rowheaderscrolled(self, first, last):
	if self._textboxLastY != first:
	    self._textboxLastY = first
	    self._textbox.yview('moveto', first)

    def _toggleHorizScrollbar(self):

	self._horizScrollbarOn = not self._horizScrollbarOn

	interior = self.interior()
	if self._horizScrollbarOn:
	    self._horizScrollbar.grid(row = 6, column = 4, sticky = 'news')
	    interior.grid_rowconfigure(5, minsize = self['scrollmargin'])
	else:
	    self._horizScrollbar.grid_forget()
	    interior.grid_rowconfigure(5, minsize = 0)

    def _toggleVertScrollbar(self):

	self._vertScrollbarOn = not self._vertScrollbarOn

	interior = self.interior()
	if self._vertScrollbarOn:
	    self._vertScrollbar.grid(row = 4, column = 6, sticky = 'news')
	    interior.grid_columnconfigure(5, minsize = self['scrollmargin'])
	else:
	    self._vertScrollbar.grid_forget()
	    interior.grid_columnconfigure(5, minsize = 0)

    # Need to explicitly forward this to override the stupid
    # (grid_)bbox method inherited from Tkinter.Frame.Grid.
    def bbox(self, index):
	return self._textbox.bbox(index)

forwardmethods(ScrolledText, Tkinter.Text, '_textbox')

######################################################################
### File: PmwHistoryText.py


_ORIGINAL = 0
_MODIFIED = 1
_DISPLAY = 2

class HistoryText(ScrolledText):

    def __init__(self, parent = None, **kw):

        # Define the megawidget options.
        optiondefs = (
	    ('compressany',         1,          None),
	    ('compresstail',        1,          None),
            ('historycommand',      None,       None),
        )
        self.defineoptions(kw, optiondefs)

        # Initialise the base class (after defining the options).
        ScrolledText.__init__(self, parent)

        # Initialise instance variables.
	self._list = []
	self._currIndex = 0
	self._pastIndex = None
	self._lastIndex = 0          # pointer to end of history list

        # Check keywords and initialise options.
        self.initialiseoptions()

    def addhistory(self):
	text = self.get()
	if text[-1] == '\n':
	    text = text[:-1]

	if len(self._list) == 0:
            # This is the first history entry.  Add it.
	    self._list.append([text, text, _MODIFIED])
            return

        currentEntry =  self._list[self._currIndex]
        if text == currentEntry[_ORIGINAL]:
            # The current history entry has not been modified. Check if
            # we need to add it again.

            if self['compresstail'] and self._currIndex == self._lastIndex:
                return

            if self['compressany']:
                return

        # Undo any changes for the current history entry, since they
        # will now be available in the new entry.
        currentEntry[_MODIFIED] = currentEntry[_ORIGINAL]

        historycommand = self['historycommand']
        if self._currIndex == self._lastIndex:
            # The last history entry is currently being displayed,
            # so disable the special meaning of the 'Next' button.
            self._pastIndex = None
            nextState = 'disabled'
        else:
            # A previous history entry is currently being displayed,
            # so allow the 'Next' button to go to the entry after this one.
            self._pastIndex = self._currIndex
            nextState = 'normal'
        if callable(historycommand):
            historycommand('normal', nextState)

        # Create the new history entry.
        self._list.append([text, text, _MODIFIED])

        # Move the pointer into the history entry list to the end.
        self._lastIndex = self._lastIndex + 1
        self._currIndex = self._lastIndex

    def next(self):
	if self._currIndex == self._lastIndex and self._pastIndex is None:
	    self.bell()
        else:
            self._modifyDisplay('next')

    def prev(self):
        self._pastIndex = None
	if self._currIndex == 0:
	    self.bell()
        else:
            self._modifyDisplay('prev')

    def undo(self):
	if len(self._list) != 0:
            self._modifyDisplay('undo')

    def redo(self):
	if len(self._list) != 0:
            self._modifyDisplay('redo')

    def gethistory(self):
        return self._list

    def _modifyDisplay(self, command):
        # Modify the display to show either the next or previous
        # history entry (next, prev) or the original or modified
        # version of the current history entry (undo, redo).

        # Save the currently displayed text.
        currentText = self.get()
        if currentText[-1] == '\n':
            currentText = currentText[:-1]

        currentEntry =  self._list[self._currIndex]
        if currentEntry[_DISPLAY] == _MODIFIED:
            currentEntry[_MODIFIED] = currentText
        elif currentEntry[_ORIGINAL] != currentText:
            currentEntry[_MODIFIED] = currentText
            if command in ('next', 'prev'):
                currentEntry[_DISPLAY] = _MODIFIED

        if command in ('next', 'prev'):
            prevstate = 'normal'
            nextstate = 'normal'
            if command == 'next':
                if self._pastIndex is not None:
                    self._currIndex = self._pastIndex
                    self._pastIndex = None
                self._currIndex = self._currIndex + 1
                if self._currIndex == self._lastIndex:
                    nextstate = 'disabled'
            elif command == 'prev':
                self._currIndex = self._currIndex - 1
                if self._currIndex == 0:
                    prevstate = 'disabled'
            historycommand = self['historycommand']
            if callable(historycommand):
                historycommand(prevstate, nextstate)
            currentEntry =  self._list[self._currIndex]
        else:
            if command == 'undo':
                currentEntry[_DISPLAY] = _ORIGINAL
            elif command == 'redo':
                currentEntry[_DISPLAY] = _MODIFIED

        # Display the new text.
        self.delete('1.0', 'end')
        self.insert('end', currentEntry[currentEntry[_DISPLAY]])

######################################################################
### File: PmwSelectionDialog.py
# Not Based on iwidgets version.



class SelectionDialog(Dialog):
    # Dialog window with selection list.
    
    # Dialog window displaying a list and requesting the user to
    # select one.

    def __init__(self, parent = None, **kw):
	# Define the megawidget options.
	
	optiondefs = (
	    ('borderx',     10,    INITOPT),
	    ('bordery',     10,    INITOPT),
	)
	self.defineoptions(kw, optiondefs)

	# Initialise the base class (after defining the options).
	Dialog.__init__(self, parent)

	# Create the components.
	interior = self.interior()
	aliases = (
	    ('listbox', 'scrolledlist_listbox'),
	    ('label', 'scrolledlist_label'),
	)
	self._list = self.createcomponent('scrolledlist',
		aliases, None,
		ScrolledListBox, (interior,),
		dblclickcommand = self.invoke)
	self._list.pack(side='top', expand='true', fill='both',
		padx = self['borderx'], pady = self['bordery'])

        if not kw.has_key('activatecommand'):
            # Whenever this dialog is activated, set the focus to the
            # ScrolledListBox's listbox widget.
            listbox = self.component('listbox')
            self.configure(activatecommand = listbox.focus_set)

	# Check keywords and initialise options.
	self.initialiseoptions()

    # Need to explicitly forward this to override the stupid
    # (grid_)size method inherited from Tkinter.Toplevel.Grid.
    def size(self):
	return self.component('listbox').size()

    # Need to explicitly forward this to override the stupid
    # (grid_)bbox method inherited from Tkinter.Toplevel.Grid.
    def bbox(self, index):
	return self.component('listbox').size(index)

forwardmethods(SelectionDialog, ScrolledListBox, '_list')

######################################################################
### File: PmwTextDialog.py
# A Dialog with a ScrolledText widget.



class TextDialog(Dialog):
    def __init__(self, parent = None, **kw):
	# Define the megawidget options.
	
	optiondefs = (
	    ('borderx',     10,    INITOPT),
	    ('bordery',     10,    INITOPT),
	)
	self.defineoptions(kw, optiondefs)

	# Initialise the base class (after defining the options).
	Dialog.__init__(self, parent)

	# Create the components.
	interior = self.interior()
	aliases = (
	    ('text', 'scrolledtext_text'),
	    ('label', 'scrolledtext_label'),
	)
	self._text = self.createcomponent('scrolledtext',
		aliases, None,
		ScrolledText, (interior,))
	self._text.pack(side='top', expand=1, fill='both',
		padx = self['borderx'], pady = self['bordery'])

	# Check keywords and initialise options.
	self.initialiseoptions()

    # Need to explicitly forward this to override the stupid
    # (grid_)bbox method inherited from Tkinter.Toplevel.Grid.
    def bbox(self, index):
	return self._text.bbox(index)

forwardmethods(TextDialog, ScrolledText, '_text')

######################################################################
### File: PmwTimeCounter.py
# Authors: Joe VanAndel and Greg McFarlane

import string
import sys
import time
import Tkinter


class TimeCounter(MegaWidget):
    """Up-down counter

    A TimeCounter is a single-line entry widget with Up and Down arrows
    which increment and decrement the Time value in the entry.  
    """

    def __init__(self, parent = None, **kw):

	# Define the megawidget options.
	
	optiondefs = (
	    ('autorepeat',    1,    None),
	    ('buttonaspect',  1.0,  INITOPT),
	    ('command',       None, None),
	    ('initwait',      300,  None),
	    ('labelmargin',   0,    INITOPT),
	    ('labelpos',      None, INITOPT),
	    ('max',           None, self._max),
	    ('min',           None, self._min),
	    ('padx',          0,    INITOPT),
	    ('pady',          0,    INITOPT),
	    ('repeatrate',    50,   None),
	    ('value',         None, INITOPT),
	)
	self.defineoptions(kw, optiondefs)

	# Initialise the base class (after defining the options).
	MegaWidget.__init__(self, parent)

    	self.arrowDirection = {}
	self._flag = 'stopped'
	self._timerId = None

	self._createComponents(kw)

	value = self['value']
	if value is None:
	    now = time.time()
	    value = time.strftime('%H:%M:%S', time.localtime(now))
    	self.setvalue(value)

	# Check keywords and initialise options.
	self.initialiseoptions()

    def _createComponents(self, kw):

	# Create the components.
	interior = self.interior()

	# If there is no label, put the arrows and the entry directly
	# into the interior, otherwise create a frame for them.  In
	# either case the border around the arrows and the entry will
	# be raised (but not around the label).
	if self['labelpos'] is None:
	    frame = interior
            if not kw.has_key('hull_relief'):
                frame.configure(relief = 'raised')
            if not kw.has_key('hull_borderwidth'):
                frame.configure(borderwidth = 1)
	else:
	    frame = self.createcomponent('frame',
		    (), None,
		    Tkinter.Frame, (interior,),
                    relief = 'raised', borderwidth = 1)
	    frame.grid(column=2, row=2, sticky='nsew')
	    interior.grid_columnconfigure(2, weight=1)
	    interior.grid_rowconfigure(2, weight=1)

	# Create the down arrow buttons.

	# Create the hour down arrow.
	self._downHourArrowBtn = self.createcomponent('downhourarrow',
		(), 'Arrow',
		Tkinter.Canvas, (frame,),
		width = 16, height = 16, relief = 'raised', borderwidth = 2)
    	self.arrowDirection[self._downHourArrowBtn] = 'down'
	self._downHourArrowBtn.grid(column = 0, row = 2)

	# Create the minute down arrow.
	self._downMinuteArrowBtn = self.createcomponent('downminutearrow',
		(), 'Arrow',
		Tkinter.Canvas, (frame,),
		width = 16, height = 16, relief = 'raised', borderwidth = 2)
    	self.arrowDirection[self._downMinuteArrowBtn] = 'down'
	self._downMinuteArrowBtn.grid(column = 1, row = 2)

	# Create the second down arrow.
	self._downSecondArrowBtn = self.createcomponent('downsecondarrow',
		(), 'Arrow',
		Tkinter.Canvas, (frame,),
		width = 16, height = 16, relief = 'raised', borderwidth = 2)
    	self.arrowDirection[self._downSecondArrowBtn] = 'down'
	self._downSecondArrowBtn.grid(column = 2, row = 2)

	# Create the entry fields.

	# Create the hour entry field.
	self._hourCounterEntry = self.createcomponent('hourentryfield',
		(('hourentry', 'hourentryfield_entry'),), None,
		EntryField, (frame,), validate='integer', entry_width = 2)
	self._hourCounterEntry.grid(column = 0, row = 1, sticky = 'news')

	# Create the minute entry field.
	self._minuteCounterEntry = self.createcomponent('minuteentryfield',
		(('minuteentry', 'minuteentryfield_entry'),), None,
		EntryField, (frame,), validate='integer', entry_width = 2)
	self._minuteCounterEntry.grid(column = 1, row = 1, sticky = 'news')

	# Create the second entry field.
	self._secondCounterEntry = self.createcomponent('secondentryfield',
		(('secondentry', 'secondentryfield_entry'),), None,
		EntryField, (frame,), validate='integer', entry_width = 2)
	self._secondCounterEntry.grid(column = 2, row = 1, sticky = 'news')

	# Create the up arrow buttons.

	# Create the hour up arrow.
	self._upHourArrowBtn = self.createcomponent('uphourarrow',
		(), 'Arrow',
		Tkinter.Canvas, (frame,),
		width = 16, height = 16, relief = 'raised', borderwidth = 2)
    	self.arrowDirection[self._upHourArrowBtn] = 'up'
	self._upHourArrowBtn.grid(column = 0, row = 0)

	# Create the minute up arrow.
	self._upMinuteArrowBtn = self.createcomponent('upminutearrow',
		(), 'Arrow',
		Tkinter.Canvas, (frame,),
		width = 16, height = 16, relief = 'raised', borderwidth = 2)
    	self.arrowDirection[self._upMinuteArrowBtn] = 'up'
	self._upMinuteArrowBtn.grid(column = 1, row = 0)

	# Create the second up arrow.
	self._upSecondArrowBtn = self.createcomponent('upsecondarrow',
		(), 'Arrow',
		Tkinter.Canvas, (frame,),
		width = 16, height = 16, relief = 'raised', borderwidth = 2)
    	self.arrowDirection[self._upSecondArrowBtn] = 'up'
	self._upSecondArrowBtn.grid(column = 2, row = 0)

	# Make it resize nicely.
	padx = self['padx']
	pady = self['pady']
	for col in range(3):
	    frame.grid_columnconfigure(col, weight = 1, pad = padx)
	frame.grid_rowconfigure(0, pad = pady)
	frame.grid_rowconfigure(2, pad = pady)

	frame.grid_rowconfigure(1, weight = 1)

	# Create the label.
	self.createlabel(interior)

	# Set bindings.

	# Up hour
	self._upHourArrowBtn.bind('<Configure>', 
		lambda  event, s=self,button=self._upHourArrowBtn: 
		s._drawArrow(button, 'up'))

	self._upHourArrowBtn.bind('<1>', 
    	    	lambda event, s=self,button=self._upHourArrowBtn: 
		s._countUp(button, 3600))

	self._upHourArrowBtn.bind('<Any-ButtonRelease-1>', 
		lambda event, s=self, button=self._upHourArrowBtn:
		s._stopUpDown(button))

	# Up minute
	self._upMinuteArrowBtn.bind('<Configure>', 
		lambda  event, s=self,button=self._upMinuteArrowBtn: 
		s._drawArrow(button, 'up'))
	    

	self._upMinuteArrowBtn.bind('<1>', 
    	    	lambda event, s=self,button=self._upMinuteArrowBtn: 
		s._countUp(button, 60))

	self._upMinuteArrowBtn.bind('<Any-ButtonRelease-1>', 
		lambda event, s=self, button=self._upMinuteArrowBtn:
		s._stopUpDown(button))

	# Up second
	self._upSecondArrowBtn.bind('<Configure>', 
		lambda  event, s=self,button=self._upSecondArrowBtn: 
		s._drawArrow(button, 'up'))
	    

	self._upSecondArrowBtn.bind('<1>', 
    	    	lambda event, s=self,button=self._upSecondArrowBtn: 
		s._countUp(button, 1))

	self._upSecondArrowBtn.bind('<Any-ButtonRelease-1>', 
		lambda event, s=self, button=self._upSecondArrowBtn:
		s._stopUpDown(button))

	# Down hour
	self._downHourArrowBtn.bind('<Configure>', 
		lambda  event, s=self,button=self._downHourArrowBtn: 
		s._drawArrow(button, 'down'))

	self._downHourArrowBtn.bind('<1>', 
    	    	lambda event, s=self,button=self._downHourArrowBtn: 
		s._countDown(button, 3600))
	self._downHourArrowBtn.bind('<Any-ButtonRelease-1>', 
		lambda event, s=self, button=self._downHourArrowBtn:
		s._stopUpDown(button))


	# Down minute
	self._downMinuteArrowBtn.bind('<Configure>', 
		lambda  event, s=self,button=self._downMinuteArrowBtn: 
		s._drawArrow(button, 'down'))

	self._downMinuteArrowBtn.bind('<1>', 
    	    	lambda event, s=self,button=self._downMinuteArrowBtn:
                s._countDown(button, 60))
	self._downMinuteArrowBtn.bind('<Any-ButtonRelease-1>', 
		lambda event, s=self, button=self._downMinuteArrowBtn:
		s._stopUpDown(button))

	# Down second
	self._downSecondArrowBtn.bind('<Configure>', 
		lambda  event, s=self,button=self._downSecondArrowBtn: 
		s._drawArrow(button, 'down'))

	self._downSecondArrowBtn.bind('<1>', 
    	    	lambda event, s=self, button=self._downSecondArrowBtn: 
		s._countDown(button,1))
	self._downSecondArrowBtn.bind('<Any-ButtonRelease-1>', 
		lambda event, s=self, button=self._downSecondArrowBtn:
		s._stopUpDown(button))

	self._hourCounterEntry.component('entry').bind(
                '<Return>', self._invoke)
	self._minuteCounterEntry.component('entry').bind(
        	'<Return>', self._invoke)
	self._secondCounterEntry.component('entry').bind(
        	'<Return>', self._invoke)

	self._hourCounterEntry.bind('<Configure>', self._resizeArrow)
	self._minuteCounterEntry.bind('<Configure>', self._resizeArrow)
	self._secondCounterEntry.bind('<Configure>', self._resizeArrow)

    def _drawArrow(self, arrow, direction):
        drawarrow(arrow, self['hourentry_foreground'], direction, 'arrow')

    def _resizeArrow(self, event = None):
	for btn in (self._upHourArrowBtn, self._upMinuteArrowBtn,
		self._upSecondArrowBtn,
		self._downHourArrowBtn,
		self._downMinuteArrowBtn, self._downSecondArrowBtn):
	    bw = (string.atoi(btn['borderwidth']) +
		    string.atoi(btn['highlightthickness']))
	    newHeight = self._hourCounterEntry.winfo_reqheight() - 2 * bw
	    newWidth = int(newHeight * self['buttonaspect'])
	    btn.configure(width=newWidth, height=newHeight)
	    self._drawArrow(btn, self.arrowDirection[btn])

    def _min(self):
	min = self['min']
        if min is None:
	    self._minVal = 0
	else:
	    self._minVal = timestringtoseconds(min)

    def _max(self):
	max = self['max']
	if max is None:
	    self._maxVal = None
	else:
	    self._maxVal = timestringtoseconds(max)

    def getvalue(self):
        return self.getstring()

    def setvalue(self, text):
        list = string.split(text, ':')
	if len(list) != 3:
	    raise ValueError('invalid value: ' + text)

	self._hour = string.atoi(list[0])
	self._minute = string.atoi(list[1])
	self._second = string.atoi(list[2]) 

    	self._setHMS()

    def getstring(self):
    	return '%02d:%02d:%02d' % (self._hour, self._minute, self._second)

    def getint(self):
    	return self._hour * 3600 + self._minute * 60 + self._second

    def _countUp(self, button, increment):
	self._relief = self._upHourArrowBtn.cget('relief')
	button.configure(relief='sunken')
	self._count(1, 'start', increment)

    def _countDown(self, button, increment):

	self._relief = self._downHourArrowBtn.cget('relief')
	button.configure(relief='sunken')
	self._count(-1, 'start', increment)

    def increment(self, seconds = 1):
	self._count(1, 'force', seconds)

    def decrement(self, seconds = 1):
	self._count(-1, 'force', seconds)

    def _count(self, factor, newFlag = None, increment = 1):
	if newFlag != 'force':
	  if newFlag is not None:
	    self._flag = newFlag

	  if self._flag == 'stopped':
	    return

	value = (string.atoi(self._hourCounterEntry.get()) *3600) + \
	      (string.atoi(self._minuteCounterEntry.get()) *60) + \
	      string.atoi(self._secondCounterEntry.get()) + \
	      factor * increment
	min = self._minVal
	max = self._maxVal
	if value < min:
	  value = min
	if max is not None and value > max:
	  value = max

	self._hour = value /3600
	self._minute = (value - (self._hour*3600)) / 60
	self._second = value - (self._hour*3600) - (self._minute*60)
	self._setHMS()

	if newFlag != 'force':
	  if self['autorepeat']:
	    if self._flag == 'start':
	      delay = self['initwait']
	      self._flag = 'running'
	    else:
	      delay = self['repeatrate']
	    self._timerId = self.after(
		delay, lambda self=self, factor=factor,increment=increment: 
		  self._count(factor,'running', increment))

    def _setHMS(self):
        self._hourCounterEntry.setentry('%02d' % self._hour)
        self._minuteCounterEntry.setentry('%02d' % self._minute)
        self._secondCounterEntry.setentry('%02d' % self._second)

    def _stopUpDown(self, button):
        if self._timerId is not None:
            self.after_cancel(self._timerId)
	    self._timerId = None
        button.configure(relief=self._relief)
        self._flag = 'stopped'

    def _invoke(self, event):
        cmd = self['command']
        if callable(cmd):
	    cmd()

    def invoke(self):
        cmd = self['command']
        if callable(cmd):
	    return cmd()

    def destroy(self):
        if self._timerId is not None:
            self.after_cancel(self._timerId)
	    self._timerId = None
        MegaWidget.destroy(self)

######################################################################
### File: PmwAboutDialog.py


class AboutDialog(MessageDialog):
    # Window to display version and contact information.

    # Class members containing resettable 'default' values:
    _version = ''
    _copyright = ''
    _contact = ''

    def __init__(self, parent = None, **kw):

	# Define the megawidget options.
	
	optiondefs = (
	    ('applicationname',   '',          INITOPT),
	    ('iconpos',           'w',         None),
	    ('icon_bitmap',       'info',      None),
	    ('buttons',           ('Close',),  None),
	    ('defaultbutton',     0,           None),
	)
	self.defineoptions(kw, optiondefs)

	# Initialise the base class (after defining the options).
	MessageDialog.__init__(self, parent)

	applicationname = self['applicationname']
        if not kw.has_key('title'):
            self.configure(title = 'About ' + applicationname)

        if not kw.has_key('message_text'):
            text = applicationname + '\n\n'
            if AboutDialog._version != '':
              text = text + 'Version ' + AboutDialog._version + '\n'
            if AboutDialog._copyright != '':
              text = text + AboutDialog._copyright + '\n\n'
            if AboutDialog._contact != '':
              text = text + AboutDialog._contact

            self.configure(message_text=text)

	# Check keywords and initialise options.
	self.initialiseoptions()

def aboutversion(value):
    AboutDialog._version = value

def aboutcopyright(value):
    AboutDialog._copyright = value

def aboutcontact(value):
    AboutDialog._contact = value

######################################################################
### File: PmwComboBox.py
# Based on iwidgets2.2.0/combobox.itk code.

import os
import string
import types
import Tkinter


class ComboBox(MegaWidget):
    def __init__(self, parent = None, **kw):

	# Define the megawidget options.
	
	optiondefs = (
	    ('autoclear',          0,          INITOPT),
	    ('buttonaspect',       1.0,        INITOPT),
	    ('dropdown',           1,          INITOPT),
	    ('fliparrow',          0,          INITOPT),
	    ('history',            1,          INITOPT),
	    ('labelmargin',        0,          INITOPT),
	    ('labelpos',           None,       INITOPT),
	    ('listheight',         200,        INITOPT),
	    ('selectioncommand',   None,       None),
	    ('sticky',            'ew',        INITOPT),
	    ('unique',             1,          INITOPT),
	)
	self.defineoptions(kw, optiondefs)

	# Initialise the base class (after defining the options).
	MegaWidget.__init__(self, parent)

	# Create the components.
	interior = self.interior()

	self._entryfield = self.createcomponent('entryfield',
		(('entry', 'entryfield_entry'),), None,
		EntryField, (interior,))
	self._entryfield.grid(column=2, row=2, sticky=self['sticky'])
	interior.grid_columnconfigure(2, weight = 1)
	self._entryWidget = self._entryfield.component('entry')

	if self['dropdown']:
	    self._isPosted = 0
            interior.grid_rowconfigure(2, weight = 1)

	    # Create the arrow button.
	    self._arrowBtn = self.createcomponent('arrowbutton',
		    (), None,
		    Tkinter.Canvas, (interior,), borderwidth = 2,
		    relief = 'raised',
		    width = 16, height = 16)
            if 'n' in self['sticky']:
                sticky = 'n'
            else:
                sticky = ''
            if 's' in self['sticky']:
                sticky = sticky + 's'
	    self._arrowBtn.grid(column=3, row=2, sticky = sticky)
	    self._arrowRelief = self._arrowBtn.cget('relief')

	    # Create the label.
	    self.createlabel(interior, childCols=2)

	    # Create the dropdown window.
	    self._popup = self.createcomponent('popup',
		    (), None,
		    Tkinter.Toplevel, (interior,))
	    self._popup.withdraw()
	    self._popup.overrideredirect(1)

	    # Create the scrolled listbox inside the dropdown window.
	    self._list = self.createcomponent('scrolledlist',
		    (('listbox', 'scrolledlist_listbox'),), None,
		    ScrolledListBox, (self._popup,),
		    hull_borderwidth = 2,
		    hull_relief = 'raised',
		    hull_height = self['listheight'],
		    usehullsize = 1,
		    listbox_exportselection = 0)
	    self._list.pack(expand=1, fill='both')
	    self.__listbox = self._list.component('listbox')

	    # Bind events to the arrow button.
	    self._arrowBtn.bind('<1>', self._postList)
	    self._arrowBtn.bind('<Configure>', self._drawArrow)
	    self._arrowBtn.bind('<3>', self._next)
	    self._arrowBtn.bind('<Shift-3>', self._previous)
	    self._arrowBtn.bind('<Down>', self._next)
	    self._arrowBtn.bind('<Up>', self._previous)
	    self._arrowBtn.bind('<Control-n>', self._next)
	    self._arrowBtn.bind('<Control-p>', self._previous)
	    self._arrowBtn.bind('<Shift-Down>', self._postList)
	    self._arrowBtn.bind('<Shift-Up>', self._postList)
	    self._arrowBtn.bind('<F34>', self._postList)
	    self._arrowBtn.bind('<F28>', self._postList)
	    self._arrowBtn.bind('<space>', self._postList)

	    # Bind events to the dropdown window.
	    self._popup.bind('<Escape>', self._unpostList)
	    self._popup.bind('<space>', self._selectUnpost)
	    self._popup.bind('<Return>', self._selectUnpost)
	    self._popup.bind('<ButtonRelease-1>', self._dropdownBtnRelease)
	    self._popup.bind('<ButtonPress-1>', self._unpostOnNextRelease)

	    # Bind events to the Tk listbox.
	    self.__listbox.bind('<Enter>', self._unpostOnNextRelease)

	    # Bind events to the Tk entry widget.
	    self._entryWidget.bind('<Configure>', self._resizeArrow)
	    self._entryWidget.bind('<Shift-Down>', self._postList)
	    self._entryWidget.bind('<Shift-Up>', self._postList)
	    self._entryWidget.bind('<F34>', self._postList)
	    self._entryWidget.bind('<F28>', self._postList)

            # Need to unpost the popup if the entryfield is unmapped (eg: 
            # its toplevel window is withdrawn) while the popup list is
            # displayed.
            self._entryWidget.bind('<Unmap>', self._unpostList)

	else:
	    # Create the scrolled listbox below the entry field.
	    self._list = self.createcomponent('scrolledlist',
		    (('listbox', 'scrolledlist_listbox'),), None,
		    ScrolledListBox, (interior,),
                    selectioncommand = self._selectCmd)
	    self._list.grid(column=2, row=3, sticky='nsew')
	    self.__listbox = self._list.component('listbox')

	    # The scrolled listbox should expand vertically.
	    interior.grid_rowconfigure(3, weight = 1)

	    # Create the label.
	    self.createlabel(interior, childRows=2)

	self._entryWidget.bind('<Down>', self._next)
	self._entryWidget.bind('<Up>', self._previous)
	self._entryWidget.bind('<Control-n>', self._next)
	self._entryWidget.bind('<Control-p>', self._previous)
	self.__listbox.bind('<Control-n>', self._next)
	self.__listbox.bind('<Control-p>', self._previous)

	if self['history']:
	    self._entryfield.configure(command=self._addHistory)

	# Check keywords and initialise options.
	self.initialiseoptions()

    def destroy(self):
	if self['dropdown'] and self._isPosted:
            popgrab(self._popup)
        MegaWidget.destroy(self)

    #======================================================================

    # Public methods

    def get(self, first = None, last=None):
	if first is None:
	    return self._entryWidget.get()
	else:
	    return self._list.get(first, last)

    def invoke(self):
	if self['dropdown']:
	    self._postList()
	else:
	    return self._selectCmd()

    def selectitem(self, index, setentry=1):
        if isinstance(index, basestring):
	    text = index
	    items = self._list.get(0, 'end')
	    if text in items:
		index = list(items).index(text)
	    else:
	    	raise IndexError('index "%s" not found' % text)
	elif setentry:
	    text = self._list.get(0, 'end')[index]

	self._list.select_clear(0, 'end')
	self._list.select_set(index, index)
	self._list.activate(index)
	self.see(index)
	if setentry:
	    self._entryfield.setentry(text)

    # Need to explicitly forward this to override the stupid
    # (grid_)size method inherited from Tkinter.Frame.Grid.
    def size(self):
	return self._list.size()

    # Need to explicitly forward this to override the stupid
    # (grid_)bbox method inherited from Tkinter.Frame.Grid.
    def bbox(self, index):
	return self._list.bbox(index)

    def clear(self):
	self._entryfield.clear()
	self._list.clear()

    #======================================================================

    # Private methods for both dropdown and simple comboboxes.

    def _addHistory(self):
	input = self._entryWidget.get()

	if input != '':
	    index = None
	    if self['unique']:
		# If item is already in list, select it and return.
		items = self._list.get(0, 'end')
		if input in items:
		    index = list(items).index(input)

	    if index is None:
		index = self._list.index('end')
		self._list.insert('end', input)

	    self.selectitem(index)
	    if self['autoclear']:
		self._entryWidget.delete(0, 'end')

	    # Execute the selectioncommand on the new entry.
	    self._selectCmd()

    def _next(self, event):
	size = self.size()
	if size <= 1:
	    return

	cursels = self.curselection()

	if len(cursels) == 0:
	    index = 0
	else:
	    index = string.atoi(cursels[0])
	    if index == size - 1:
		index = 0
	    else:
		index = index + 1

	self.selectitem(index)

    def _previous(self, event):
	size = self.size()
	if size <= 1:
	    return

	cursels = self.curselection()

	if len(cursels) == 0:
	    index = size - 1
	else:
	    index = string.atoi(cursels[0])
	    if index == 0:
		index = size - 1
	    else:
		index = index - 1

	self.selectitem(index)

    def _selectCmd(self, event=None):

	sels = self.getcurselection()
	if len(sels) == 0:
	    item = None
	else:
	    item = sels[0]
	    self._entryfield.setentry(item)

	cmd = self['selectioncommand']
	if callable(cmd):
            if event is None:
                # Return result of selectioncommand for invoke() method.
                return cmd(item)
            else:
                cmd(item)

    #======================================================================

    # Private methods for dropdown combobox.

    def _drawArrow(self, event=None, sunken=0):
        arrow = self._arrowBtn
	if sunken:
	    self._arrowRelief = arrow.cget('relief')
	    arrow.configure(relief = 'sunken')
	else:
	    arrow.configure(relief = self._arrowRelief)

	if self._isPosted and self['fliparrow']:
            direction = 'up'
        else:
            direction = 'down'
        drawarrow(arrow, self['entry_foreground'], direction, 'arrow')

    def _postList(self, event = None):
        self._isPosted = 1
        self._drawArrow(sunken=1)

        # Make sure that the arrow is displayed sunken.
        self.update_idletasks()

        x = self._entryfield.winfo_rootx()
        y = self._entryfield.winfo_rooty() + \
            self._entryfield.winfo_height()
        w = self._entryfield.winfo_width() + self._arrowBtn.winfo_width()
        h =  self.__listbox.winfo_height()
        sh = self.winfo_screenheight()

        if y + h > sh and y > sh / 2:
            y = self._entryfield.winfo_rooty() - h

        self._list.configure(hull_width=w)

        setgeometryanddeiconify(self._popup, '+%d+%d' % (x, y))

        # Grab the popup, so that all events are delivered to it, and
        # set focus to the listbox, to make keyboard navigation
        # easier.
        pushgrab(self._popup, 1, self._unpostList)
        self.__listbox.focus_set()

        self._drawArrow()

        # Ignore the first release of the mouse button after posting the
        # dropdown list, unless the mouse enters the dropdown list.
        self._ignoreRelease = 1

    def _dropdownBtnRelease(self, event):
	if (event.widget == self._list.component('vertscrollbar') or
		event.widget == self._list.component('horizscrollbar')):
	    return

	if self._ignoreRelease:
	    self._unpostOnNextRelease()
	    return

        self._unpostList()

	if (event.x >= 0 and event.x < self.__listbox.winfo_width() and
		event.y >= 0 and event.y < self.__listbox.winfo_height()):
	    self._selectCmd()

    def _unpostOnNextRelease(self, event = None):
	self._ignoreRelease = 0

    def _resizeArrow(self, event):
	bw = (string.atoi(self._arrowBtn['borderwidth']) + 
		string.atoi(self._arrowBtn['highlightthickness']))
	newHeight = self._entryfield.winfo_reqheight() - 2 * bw
	newWidth = int(newHeight * self['buttonaspect'])
	self._arrowBtn.configure(width=newWidth, height=newHeight)
	self._drawArrow()

    def _unpostList(self, event=None):
	if not self._isPosted:
            # It is possible to get events on an unposted popup.  For
            # example, by repeatedly pressing the space key to post
            # and unpost the popup.  The <space> event may be
            # delivered to the popup window even though
            # popgrab() has set the focus away from the
            # popup window.  (Bug in Tk?)
            return

        # Restore the focus before withdrawing the window, since
        # otherwise the window manager may take the focus away so we
        # can't redirect it.  Also, return the grab to the next active
        # window in the stack, if any.
        popgrab(self._popup)
	self._popup.withdraw()

	self._isPosted = 0
	self._drawArrow()

    def _selectUnpost(self, event):
        self._unpostList()
	self._selectCmd()

forwardmethods(ComboBox, ScrolledListBox, '_list')
forwardmethods(ComboBox, EntryField, '_entryfield')

######################################################################
### File: PmwComboBoxDialog.py
# Not Based on iwidgets version.



class ComboBoxDialog(Dialog):
    # Dialog window with simple combobox.
    
    # Dialog window displaying a list and entry field and requesting
    # the user to make a selection or enter a value

    def __init__(self, parent = None, **kw):
	# Define the megawidget options.
	
	optiondefs = (
	    ('borderx',    10,              INITOPT),
	    ('bordery',    10,              INITOPT),
	)
	self.defineoptions(kw, optiondefs)

	# Initialise the base class (after defining the options).
	Dialog.__init__(self, parent)

	# Create the components.
	interior = self.interior()

	aliases = (
	    ('listbox', 'combobox_listbox'),
	    ('scrolledlist', 'combobox_scrolledlist'),
	    ('entry', 'combobox_entry'),
	    ('label', 'combobox_label'),
	)
	self._combobox = self.createcomponent('combobox',
		aliases, None,
		ComboBox, (interior,),
		scrolledlist_dblclickcommand = self.invoke,
		dropdown = 0,
	)
	self._combobox.pack(side='top', expand='true', fill='both',
		padx = self['borderx'], pady = self['bordery'])

        if not kw.has_key('activatecommand'):
            # Whenever this dialog is activated, set the focus to the
            # ComboBox's listbox widget.
            listbox = self.component('listbox')
            self.configure(activatecommand = listbox.focus_set)

	# Check keywords and initialise options.
	self.initialiseoptions()

    # Need to explicitly forward this to override the stupid
    # (grid_)size method inherited from Tkinter.Toplevel.Grid.
    def size(self):
	return self._combobox.size()

    # Need to explicitly forward this to override the stupid
    # (grid_)bbox method inherited from Tkinter.Toplevel.Grid.
    def bbox(self, index):
	return self._combobox.bbox(index)

forwardmethods(ComboBoxDialog, ComboBox, '_combobox')

######################################################################
### File: PmwCounter.py
import string
import sys
import types
import Tkinter


class Counter(MegaWidget):

    def __init__(self, parent = None, **kw):

	# Define the megawidget options.
	
	optiondefs = (
	    ('autorepeat',     1,             None),
	    ('buttonaspect',   1.0,           INITOPT),
	    ('datatype',       'numeric',     self._datatype),
	    ('increment',      1,             None),
	    ('initwait',       300,           None),
	    ('labelmargin',    0,             INITOPT),
	    ('labelpos',       None,          INITOPT),
	    ('orient',         'horizontal',  INITOPT),
	    ('padx',           0,             INITOPT),
	    ('pady',           0,             INITOPT),
	    ('repeatrate',     50,            None),
	    ('sticky',         'ew',          INITOPT),
	)
	self.defineoptions(kw, optiondefs)

	# Initialise the base class (after defining the options).
	MegaWidget.__init__(self, parent)

	# Initialise instance variables.
	self._timerId = None
	self._normalRelief = None

	# Create the components.
	interior = self.interior()

	# If there is no label, put the arrows and the entry directly
	# into the interior, otherwise create a frame for them.  In
	# either case the border around the arrows and the entry will
	# be raised (but not around the label).
	if self['labelpos'] is None:
	    frame = interior
            if not kw.has_key('hull_relief'):
                frame.configure(relief = 'raised')
            if not kw.has_key('hull_borderwidth'):
                frame.configure(borderwidth = 1)
	else:
	    frame = self.createcomponent('frame',
		    (), None,
		    Tkinter.Frame, (interior,),
                    relief = 'raised', borderwidth = 1)
	    frame.grid(column=2, row=2, sticky=self['sticky'])
	    interior.grid_columnconfigure(2, weight=1)
	    interior.grid_rowconfigure(2, weight=1)

	# Create the down arrow.
	self._downArrowBtn = self.createcomponent('downarrow',
		(), 'Arrow',
		Tkinter.Canvas, (frame,),
		width = 16, height = 16, relief = 'raised', borderwidth = 2)

	# Create the entry field.
	self._counterEntry = self.createcomponent('entryfield',
		(('entry', 'entryfield_entry'),), None,
		EntryField, (frame,))

	# Create the up arrow.
	self._upArrowBtn = self.createcomponent('uparrow',
		(), 'Arrow',
		Tkinter.Canvas, (frame,),
		width = 16, height = 16, relief = 'raised', borderwidth = 2)

	padx = self['padx']
	pady = self['pady']
	orient = self['orient']
	if orient == 'horizontal':
	    self._downArrowBtn.grid(column = 0, row = 0)
	    self._counterEntry.grid(column = 1, row = 0,
                    sticky = self['sticky'])
	    self._upArrowBtn.grid(column = 2, row = 0)
	    frame.grid_columnconfigure(1, weight = 1)
	    frame.grid_rowconfigure(0, weight = 1)
	    if Tkinter.TkVersion >= 4.2:
		frame.grid_columnconfigure(0, pad = padx)
		frame.grid_columnconfigure(2, pad = padx)
		frame.grid_rowconfigure(0, pad = pady)
	elif orient == 'vertical':
	    self._upArrowBtn.grid(column = 0, row = 0, sticky = 's')
	    self._counterEntry.grid(column = 0, row = 1,
                    sticky = self['sticky'])
	    self._downArrowBtn.grid(column = 0, row = 2, sticky = 'n')
	    frame.grid_columnconfigure(0, weight = 1)
	    frame.grid_rowconfigure(0, weight = 1)
	    frame.grid_rowconfigure(2, weight = 1)
	    if Tkinter.TkVersion >= 4.2:
		frame.grid_rowconfigure(0, pad = pady)
		frame.grid_rowconfigure(2, pad = pady)
		frame.grid_columnconfigure(0, pad = padx)
	else:
	    raise ValueError('bad orient option ' + repr(orient)
                             + ': must be either \'horizontal\' or \'vertical\'')

	self.createlabel(interior)

	self._upArrowBtn.bind('<Configure>', self._drawUpArrow)
	self._upArrowBtn.bind('<1>', self._countUp)
	self._upArrowBtn.bind('<Any-ButtonRelease-1>', self._stopCounting)
	self._downArrowBtn.bind('<Configure>', self._drawDownArrow)
	self._downArrowBtn.bind('<1>', self._countDown)
	self._downArrowBtn.bind('<Any-ButtonRelease-1>', self._stopCounting)
	self._counterEntry.bind('<Configure>', self._resizeArrow)
	entry = self._counterEntry.component('entry')
	entry.bind('<Down>', lambda event, s = self: s._key_decrement(event))
	entry.bind('<Up>', lambda event, s = self: s._key_increment(event))

	# Need to cancel the timer if an arrow button is unmapped (eg: 
	# its toplevel window is withdrawn) while the mouse button is
	# held down.  The canvas will not get the ButtonRelease event
	# if it is not mapped, since the implicit grab is cancelled.
	self._upArrowBtn.bind('<Unmap>', self._stopCounting)
	self._downArrowBtn.bind('<Unmap>', self._stopCounting)

	# Check keywords and initialise options.
	self.initialiseoptions()

    def _resizeArrow(self, event):
	for btn in (self._upArrowBtn, self._downArrowBtn):
	    bw = (string.atoi(btn['borderwidth']) +
		    string.atoi(btn['highlightthickness']))
	    newHeight = self._counterEntry.winfo_reqheight() - 2 * bw
	    newWidth = int(newHeight * self['buttonaspect'])
	    btn.configure(width=newWidth, height=newHeight)
	    self._drawArrow(btn)

    def _drawUpArrow(self, event):
	self._drawArrow(self._upArrowBtn)

    def _drawDownArrow(self, event):
	self._drawArrow(self._downArrowBtn)

    def _drawArrow(self, arrow):
        if self['orient'] == 'vertical':
            if arrow == self._upArrowBtn:
                direction = 'up'
            else:
                direction = 'down'
        else:
            if arrow == self._upArrowBtn:
                direction = 'right'
            else:
                direction = 'left'
        drawarrow(arrow, self['entry_foreground'], direction, 'arrow')

    def _stopCounting(self, event = None):
        if self._timerId is not None:
            self.after_cancel(self._timerId)
	    self._timerId = None
	if self._normalRelief is not None:
	    button, relief = self._normalRelief
	    button.configure(relief=relief)
	    self._normalRelief = None

    def _countUp(self, event):
	self._normalRelief = (self._upArrowBtn, self._upArrowBtn.cget('relief'))
	self._upArrowBtn.configure(relief='sunken')
	# Force arrow down (it may come up immediately, if increment fails).
	self._upArrowBtn.update_idletasks()
	self._count(1, 1)

    def _countDown(self, event):
	self._normalRelief = (self._downArrowBtn, self._downArrowBtn.cget('relief'))
	self._downArrowBtn.configure(relief='sunken')
	# Force arrow down (it may come up immediately, if increment fails).
	self._downArrowBtn.update_idletasks()
	self._count(-1, 1)

    def increment(self):
	self._forceCount(1)

    def decrement(self):
	self._forceCount(-1)

    def _key_increment(self, event):
	self._forceCount(1)
	self.update_idletasks()

    def _key_decrement(self, event):
	self._forceCount(-1)
	self.update_idletasks()

    def _datatype(self):
	datatype = self['datatype']

        if isinstance(datatype, dict):
	    self._counterArgs = datatype.copy()
	    if self._counterArgs.has_key('counter'):
		datatype = self._counterArgs['counter']
		del self._counterArgs['counter']
	    else:
		datatype = 'numeric'
	else:
	    self._counterArgs = {}

	if _counterCommands.has_key(datatype):
	    self._counterCommand = _counterCommands[datatype]
	elif callable(datatype):
	    self._counterCommand = datatype
	else:
	    validValues = _counterCommands.keys()
	    validValues.sort()
	    raise ValueError('bad datatype value "%s":  must be a function '
                             'or one of %s' % (datatype, validValues))

    def _forceCount(self, factor):
	if not self.valid():
	    self.bell()
	    return

	text = self._counterEntry.get()
	try:
	    value = apply(self._counterCommand,
		    (text, factor, self['increment']), self._counterArgs)
	except ValueError:
	    self.bell()
	    return

        previousICursor = self._counterEntry.index('insert')
	if self._counterEntry.setentry(value) == OK:
	    self._counterEntry.xview('end')
	    self._counterEntry.icursor(previousICursor)

    def _count(self, factor, first):
	if not self.valid():
	    self.bell()
	    return

	self._timerId = None
	origtext = self._counterEntry.get()
	try:
	    value = apply(self._counterCommand,
		    (origtext, factor, self['increment']), self._counterArgs)
	except ValueError:
	    # If text is invalid, stop counting.
	    self._stopCounting()
	    self.bell()
	    return

	# If incrementing produces an invalid value, restore previous
	# text and stop counting.
        previousICursor = self._counterEntry.index('insert')
	valid = self._counterEntry.setentry(value)
	if valid != OK:
	    self._stopCounting()
	    self._counterEntry.setentry(origtext)
	    if valid == PARTIAL:
		self.bell()
	    return
	self._counterEntry.xview('end')
	self._counterEntry.icursor(previousICursor)

	if self['autorepeat']:
	    if first:
		delay = self['initwait']
	    else:
		delay = self['repeatrate']
	    self._timerId = self.after(delay,
		    lambda self=self, factor=factor: self._count(factor, 0))

    def destroy(self):
	self._stopCounting()
        MegaWidget.destroy(self)

forwardmethods(Counter, EntryField, '_counterEntry')

def _changeNumber(text, factor, increment):
  value = string.atol(text)
  if factor > 0:
    value = (value / increment) * increment + increment
  else:
    value = ((value - 1) / increment) * increment

  # Get rid of the 'L' at the end of longs (in python up to 1.5.2).
  rtn = str(value)
  if rtn[-1] == 'L':
      return rtn[:-1]
  else:
      return rtn

def _changeReal(text, factor, increment, separator = '.'):
  value = stringtoreal(text, separator)
  div = value / increment

  # Compare reals using str() to avoid problems caused by binary
  # numbers being only approximations to decimal numbers.
  # For example, if value is -0.3 and increment is 0.1, then
  # int(value/increment) = -2, not -3 as one would expect.
  if str(div)[-2:] == '.0':
    # value is an even multiple of increment.
    div = round(div) + factor
  else:
    # value is not an even multiple of increment.
    div = int(div) * 1.0
    if value < 0:
      div = div - 1
    if factor > 0:
      div = (div + 1)

  value = div * increment

  text = str(value)
  if separator != '.':
      index = string.find(text, '.')
      if index >= 0:
	text = text[:index] + separator + text[index + 1:]
  return text

def _changeDate(value, factor, increment, format = 'ymd',
	separator = '/', yyyy = 0):

  jdn = datestringtojdn(value, format, separator) + factor * increment

  y, m, d = jdntoymd(jdn)
  result = ''
  for index in range(3):
    if index > 0:
      result = result + separator
    f = format[index]
    if f == 'y':
      if yyyy:
        result = result + '%02d' % y
      else:
        result = result + '%02d' % (y % 100)
    elif f == 'm':
      result = result + '%02d' % m
    elif f == 'd':
      result = result + '%02d' % d

  return result

_SECSPERDAY = 24 * 60 * 60
def _changeTime(value, factor, increment, separator = ':', time24 = 0):
  unixTime = timestringtoseconds(value, separator)
  if factor > 0:
    chunks = unixTime / increment + 1
  else:
    chunks = (unixTime - 1) / increment
  unixTime = chunks * increment
  if time24:
      while unixTime < 0:
	  unixTime = unixTime + _SECSPERDAY
      while unixTime >= _SECSPERDAY:
	  unixTime = unixTime - _SECSPERDAY
  if unixTime < 0:
    unixTime = -unixTime
    sign = '-'
  else:
    sign = ''
  secs = unixTime % 60
  unixTime = unixTime / 60
  mins = unixTime % 60
  hours = unixTime / 60
  return '%s%02d%s%02d%s%02d' % (sign, hours, separator, mins, separator, secs)

# hexadecimal, alphabetic, alphanumeric not implemented
_counterCommands = {
    'numeric'   : _changeNumber,      # } integer
    'integer'   : _changeNumber,      # } these two use the same function
    'real'      : _changeReal,        # real number
    'time'      : _changeTime,
    'date'      : _changeDate,
}

######################################################################
### File: PmwCounterDialog.py


# A Dialog with a counter

class CounterDialog(Dialog):

    def __init__(self, parent = None, **kw):

	# Define the megawidget options.
	
	optiondefs = (
	    ('borderx',    20,  INITOPT),
	    ('bordery',    20,  INITOPT),
	)
	self.defineoptions(kw, optiondefs)

	# Initialise the base class (after defining the options).
	Dialog.__init__(self, parent)

	# Create the components.
	interior = self.interior()

	# Create the counter.
	aliases = (
	    ('entryfield', 'counter_entryfield'),
	    ('entry', 'counter_entryfield_entry'),
	    ('label', 'counter_label')
	)
	self._cdCounter = self.createcomponent('counter',
		aliases, None,
		Counter, (interior,))
	self._cdCounter.pack(fill='x', expand=1,
		padx = self['borderx'], pady = self['bordery'])
	
        if not kw.has_key('activatecommand'):
            # Whenever this dialog is activated, set the focus to the
            # Counter's entry widget.
            tkentry = self.component('entry')
            self.configure(activatecommand = tkentry.focus_set)

	# Check keywords and initialise options.
	self.initialiseoptions()

    # Supply aliases to some of the entry component methods.
    def insertentry(self, index, text):
	self._cdCounter.insert(index, text)

    def deleteentry(self, first, last=None):
	self._cdCounter.delete(first, last)

    def indexentry(self, index):
	return self._cdCounter.index(index)

forwardmethods(CounterDialog, Counter, '_cdCounter')

######################################################################
### File: PmwLogicalFont.py
import os
import string

def _font_initialise(root, size=None, fontScheme = None):
    global _fontSize
    if size is not None:
        _fontSize = size

    if fontScheme in ('pmw1', 'pmw2'):
        if os.name == 'posix':
            defaultFont = logicalfont('Helvetica')
            menuFont = logicalfont('Helvetica', weight='bold', slant='italic')
            scaleFont = logicalfont('Helvetica', slant='italic')
            root.option_add('*Font',            defaultFont,  'userDefault')
            root.option_add('*Menu*Font',       menuFont,     'userDefault')
            root.option_add('*Menubutton*Font', menuFont,     'userDefault')
            root.option_add('*Scale.*Font',     scaleFont,    'userDefault')

            if fontScheme == 'pmw1':
                balloonFont = logicalfont('Helvetica', -6, pixel = '12')
            else: # fontScheme == 'pmw2'
                balloonFont = logicalfont('Helvetica', -2)
            root.option_add('*Balloon.*Font', balloonFont, 'userDefault')
        else:
            defaultFont = logicalfont('Helvetica')
            root.option_add('*Font', defaultFont,  'userDefault')
    elif fontScheme == 'default':
        defaultFont = ('Helvetica', '-%d' % (_fontSize,), 'bold')
        entryFont = ('Helvetica', '-%d' % (_fontSize,))
        textFont = ('Courier', '-%d' % (_fontSize,))
        root.option_add('*Font',            defaultFont,  'userDefault')
        root.option_add('*Entry*Font',      entryFont,    'userDefault')
        root.option_add('*Text*Font',       textFont,     'userDefault')

def logicalfont(name='Helvetica', sizeIncr = 0, **kw):
  if not _fontInfo.has_key(name):
    raise ValueError('font %s does not exist' % name)

  rtn = []
  for field in _fontFields:
    if kw.has_key(field):
      logicalValue = kw[field]
    elif _fontInfo[name].has_key(field):
      logicalValue = _fontInfo[name][field]
    else:
      logicalValue = '*'

    if _propertyAliases[name].has_key((field, logicalValue)):
      realValue = _propertyAliases[name][(field, logicalValue)]
    elif _propertyAliases[name].has_key((field, None)):
      realValue = _propertyAliases[name][(field, None)]
    elif _propertyAliases[None].has_key((field, logicalValue)):
      realValue = _propertyAliases[None][(field, logicalValue)]
    elif _propertyAliases[None].has_key((field, None)):
      realValue = _propertyAliases[None][(field, None)]
    else:
      realValue = logicalValue

    if field == 'size':
      if realValue == '*':
	  realValue = _fontSize
      realValue = str((realValue + sizeIncr) * 10)

    rtn.append(realValue)

  return string.join(rtn, '-')

def logicalfontnames():
  return _fontInfo.keys()

if os.name == 'nt':
    _fontSize = 16
else:
    _fontSize = 14

_fontFields = (
  'registry', 'foundry', 'family', 'weight', 'slant', 'width', 'style',
  'pixel', 'size', 'xres', 'yres', 'spacing', 'avgwidth', 'charset', 'encoding')

# <_propertyAliases> defines other names for which property values may
# be known by.  This is required because italics in adobe-helvetica
# are specified by 'o', while other fonts use 'i'.

_propertyAliases = {}

_propertyAliases[None] = {
  ('slant', 'italic') : 'i',
  ('slant', 'normal') : 'r',
  ('weight', 'light') : 'normal',
  ('width', 'wide') : 'normal',
  ('width', 'condensed') : 'normal',
}

# <_fontInfo> describes a 'logical' font, giving the default values of
# some of its properties.

_fontInfo = {}

_fontInfo['Helvetica'] = {
  'foundry' : 'adobe',
  'family' : 'helvetica',
  'registry' : '',
  'charset' : 'iso8859',
  'encoding' : '1',
  'spacing' : 'p',
  'slant' : 'normal',
  'width' : 'normal',
  'weight' : 'normal',
}

_propertyAliases['Helvetica'] = {
  ('slant', 'italic') : 'o',
  ('weight', 'normal') : 'medium',
  ('weight', 'light') : 'medium',
}

_fontInfo['Times'] = {
  'foundry' : 'adobe',
  'family' : 'times',
  'registry' : '',
  'charset' : 'iso8859',
  'encoding' : '1',
  'spacing' : 'p',
  'slant' : 'normal',
  'width' : 'normal',
  'weight' : 'normal',
}

_propertyAliases['Times'] = {
  ('weight', 'normal') : 'medium',
  ('weight', 'light') : 'medium',
}

_fontInfo['Fixed'] = {
  'foundry' : 'misc',
  'family' : 'fixed',
  'registry' : '',
  'charset' : 'iso8859',
  'encoding' : '1',
  'spacing' : 'c',
  'slant' : 'normal',
  'width' : 'normal',
  'weight' : 'normal',
}

_propertyAliases['Fixed'] = {
  ('weight', 'normal') : 'medium',
  ('weight', 'light') : 'medium',
  ('style', None) : '',
  ('width', 'condensed') : 'semicondensed',
}

_fontInfo['Courier'] = {
  'foundry' : 'adobe',
  'family' : 'courier',
  'registry' : '',
  'charset' : 'iso8859',
  'encoding' : '1',
  'spacing' : 'm',
  'slant' : 'normal',
  'width' : 'normal',
  'weight' : 'normal',
}

_propertyAliases['Courier'] = {
  ('weight', 'normal') : 'medium',
  ('weight', 'light') : 'medium',
  ('style', None) : '',
}

_fontInfo['Typewriter'] = {
  'foundry' : 'b&h',
  'family' : 'lucidatypewriter',
  'registry' : '',
  'charset' : 'iso8859',
  'encoding' : '1',
  'spacing' : 'm',
  'slant' : 'normal',
  'width' : 'normal',
  'weight' : 'normal',
}

_propertyAliases['Typewriter'] = {
  ('weight', 'normal') : 'medium',
  ('weight', 'light') : 'medium',
}

if os.name == 'nt':
    # For some reason 'fixed' fonts on NT aren't.
    _fontInfo['Fixed'] = _fontInfo['Courier']
    _propertyAliases['Fixed'] = _propertyAliases['Courier']

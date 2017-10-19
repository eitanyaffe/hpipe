include ../makeshift.mk

# load modules
$(call _module, module1/module1.mk)
$(call _module, module2/module2.mk)

# define class point which controls two variables x,y
$(call _class,point,x y)

# define instance p1 and values x=3,y=4
$(call _class_instance,point,p1,X 4)

X=17

# set active point to be p1
$(call _class_activate,point,p1)
$(info x=$(x))
$(info y=$(y))

# global variable
G_V1?=3

# active module
$(call _active_module,$(m))

example:
	$(MAKE) m=module1 step1
	$(MAKE) m=module2 step2

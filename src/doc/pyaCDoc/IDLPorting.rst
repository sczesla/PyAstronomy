The PyA guide for porting IDL code
======================================

The "Interactive Data Language" (IDL) has been an influential tool in the
astronomical community. Some essential code is written in IDL, and sometimes
it may be desirable to port it to Python.

Some typical aspects of IDL cannot be directly translated into Python
constructs. Therefore, we provide a number of guidelines here, which can
be used in porting IDL code to Python:

 * **Case sensitivity**
    IDL is case-insensitive. We recommend to use lower case for all function and
    variable names.
 * **IDL keywords**
    IDL functions may take keyword arguments. If the purpose of the keyword is only
    to be present or not, we recommend to use a boolean keyword that defaults to False
    in Python. If the keyword argument is supposed to convey information to the function,
    we recommend to use a Python keyword that defaults to the standard value used (in IDL) when
    when the keyword is not given.
 * **Arrays**
    Python provides no built-in array support. IDL arrays are represented by numpy arrays.
    Were IDL uses an array as return value, it may not be good practice to do the same in
    Python. In this case, it remains up to the developer to choose a reasonable solution
    and document it.
 * **Functions supporting both single float and array arguments**
    We regard array support in Python as optional. In most cases performance is not an issue for
    the functions under consideration.   
      
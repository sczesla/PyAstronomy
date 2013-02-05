PyAstronomy versioning rules
==============================

In the following, the version tag of PyAstronomy is described. \
The version tag consists of three numbers and a possible extension, \
specifying the status (beta) of the version. \

The entire PyA package is characterized by a single version number.

As an example, such a version number could look like:
 - **0.1.2**
 - **0.1.2-beta**.

The version tag, thus, consists of four parts being:
 - *First number*: An integer specifying the major release number.
 - *Second number*: An integer specifying the minor release number. A change in this number \
   indicates that an API change can have occurred (but not necessarily).
 - *Third number*: An integer specifying the fix release number. A change in this number \
   indicates that no changes to the public API, possibly breaking running code, \
   have been introduced; though, the internal working of one or more submodules \
   may have changed.
 - *-beta*: Indicates that a version is currently developed on the basis of the stated version. \
   **No guarantee for nothing**.  
                   
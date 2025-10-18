from __future__ import print_function
from . import pyaErrors as PE
import difflib
import six
 

def fuzzyMatch(inkey, wordList, caseSensitive=True, n=3, cutoff=0.6, raises=False):
  """
    Find exact and approximate matches for input in word list.
    
    Uses `get_close_matches` from Python's difflib to search
    for the best matches between the input keyword and a list
    of words.
    
    Parameters
    ----------
    inkey : string
        The input keyword.
    wordList : list of strings
        List of words with possible matches for `inkey`.
    caseSensitive : boolean, optional
        If True (default), the search will be case
        sensitive.
    n : int, optional
        Number of potential matches returned by
        `get_close_matches`.
    cutoff : float, optional
        Number between 0 and 1 indicating the degree
        of similarity between `inkey` and the entries from
        `wordList`. The lower the number, the more dissimilar
        the potential matches may be. (`cutoff` parameter
        from `get_close_matches`).
    raises : boolean, optional
        If True, a `PyAValError` giving a summary of the
        failure will be raised if no exact match is found.
        The default is false.
    
    Returns
    -------
    Matches : dictionary
        If found, contains the exact match (key "em") found
        in the list (in lower case if `caseSensitive` is True)
        and a list of close matches (key "cm"), which the user
        may have meant to specify.  If an exact match is found,
        also its index in `wordList` is provided.
  """
  if not isinstance(wordList[0], six.string_types):
    raise(PE.PyAValError("Worlist does not seem to be a list of strings. First item has type: " + str(type(wordList[0])), \
                         where="fuzzyMatch", \
                         solution="Convert into string."))
  if not caseSensitive:
    key = inkey.lower()
    wordList = [x.lower() for x in wordList]
  else:
    key = inkey
  result = {"em":None, "cm":None, "index":None}
  result["cm"] = difflib.get_close_matches(key, wordList, n, cutoff)
  if len(result["cm"]) == 0:
    if raises:
      raise(PE.PyAValError("No match for '" + str(inkey) + "'", \
            solution="Check the input keyword and/or the data."))
  if key == result["cm"][0]:
    # An exact match between the input list
    result["em"] = key
    result["index"] = wordList.index(key)
  else:
    # There is no exact, but a approximate matches
    if raises:
      raise(PE.PyAValError("No exact match found for key '" + str(inkey) + "'\n" + \
                           "Did you mean any of these: " + ', '.join(result["cm"]) + "?", \
              solution=["Check the input keyword and/or the data."]))
  return result
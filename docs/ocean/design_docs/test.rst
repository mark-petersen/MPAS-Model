.. role:: raw-math(raw)
    :format: latex html

* This is a bulleted list.
* It has two items, the second
  item uses two lines.

1. This is a numbered list.
2. It has two items too.

#. This is a numbered list.
#. It has two items too.
Nested lists are possible, but be aware that they must be separated from the parent list items by blank lines:

* this is
* a list

  * with a nested list
  * and some subitems

* and here the parent list continues
Definition lists (ref) are created as follows:

term (up to a line of text)
   Definition of the term, which must be indented

   and can even consist of multiple paragraphs

next term
   Description.
Note that the term cannot have more than one line of text.

Quoted paragraphs (ref) are created by just indenting them more than the surrounding paragraphs.

Line blocks (ref) are a way of preserving line breaks:

| These lines are
| broken exactly like in
| the source file.
There are also several more special blocks available:

field lists (ref, with caveats noted in Field Lists)
option lists (ref)
quoted literal blocks (ref)
doctest blocks (ref)
Literal blocks
Literal code blocks (ref) are introduced by ending a paragraph with the special marker ::. The literal block must be indented (and, like all paragraphs, separated from the surrounding ones by blank lines):

This is a normal text paragraph. The next paragraph is a code sample::

   It is not processed in any way, except
   that the indentation is removed.

   It can span multiple lines.

This is a normal text paragraph again.
The handling of the :: marker is smart:

If it occurs as a paragraph of its own, that paragraph is completely left out of the document.
If it is preceded by whitespace, the marker is removed.
If it is preceded by non-whitespace, the marker is replaced by a single colon.
That way, the second sentence in the above example’s first paragraph would be rendered as “The next paragraph is a code sample:”.

Code highlighting can be enabled for these literal blocks on a document-wide basis using the highlight directive and on a project-wide basis using the highlight_language configuration option. The code-block directive can be used to set highlighting on a block-by-block basis. These directives are discussed later.

Doctest blocks
Doctest blocks (ref) are interactive Python sessions cut-and-pasted into docstrings. They do not require the literal blocks syntax. The doctest block must end with a blank line and should not end with with an unused prompt:

>>> 1 + 1
2
Tables
For grid tables (ref), you have to “paint” the cell grid yourself. They look like this:

+------------------------+------------+----------+----------+
| Header row, column 1   | Header 2   | Header 3 | Header 4 |
| (header rows optional) |            |          |          |
+========================+============+==========+==========+
| body row 1, column 1   | column 2   | column 3 | column 4 |
+------------------------+------------+----------+----------+
| body row 2             | ...        | ...      |          |
+------------------------+------------+----------+----------+
Simple tables (ref) are easier to write, but limited: they must contain more than one row, and the first column cells cannot contain multiple lines. They look like this:

=====  =====  =======
A      B      A and B
=====  =====  =======
False  False  False
True   False  False
False  True   False
True   True   True
=====  =====  =======
Two more syntaxes are supported: CSV tables and List tables. They use an explicit markup block. Refer to Tables for more information.

Hyperlinks
External links
Use `Link text <https://domain.invalid/>`_ for inline web links. If the link text should be the web address, you don’t need special markup at all, the parser finds links and mail addresses in ordinary text.

Important

There must be a space between the link text and the opening < for the URL.

You can also separate the link and the target definition (ref), like this:

This is a paragraph that contains `a link`_.

.. _a link: https://domain.invalid/
Internal links
Internal linking is done via a special reST role provided by Sphinx, see the section on specific markup, Cross-referencing arbitrary locations.

Sections
Section headers (ref) are created by underlining (and optionally overlining) the section title with a punctuation character, at least as long as the text:

=================
This is a heading
=================

This is an equation.:math: `A=r^2`.
So is this ..math::
a^2 + b^2 = c^2

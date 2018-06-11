#ifndef GAPPA_TOOLS_VERSION_H_
#define GAPPA_TOOLS_VERSION_H_

/*
    gappa - Genesis Applications for Phylogenetic Placement Analysis
    Copyright (C) 2017-2018 Lucas Czech and HITS gGmbH

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

    Contact:
    Lucas Czech <lucas.czech@h-its.org>
    Exelixis Lab, Heidelberg Institute for Theoretical Studies
    Schloss-Wolfsbrunnenweg 35, D-69118 Heidelberg, Germany
*/

#include <string>

// =================================================================================================
//      Gappa Version
// =================================================================================================

inline std::string gappa_version()
{
    return "v0.0.0"; // #GAPPA_VERSION#
}

inline std::string gappa_header()
{
    return "\
                                              ....      ....  \n\
                                             '' '||.   .||'   \n\
                                                  ||  ||      \n\
                                                  '|.|'       \n\
     ...'   ....   ... ...  ... ...   ....        .|'|.       \n\
    |  ||  '' .||   ||'  ||  ||'  || '' .||      .|'  ||      \n\
     |''   .|' ||   ||    |  ||    | .|' ||     .|'|.  ||     \n\
    '....  '|..'|'. ||...'   ||...'  '|..'|.    '||'    ||:.  \n\
    '....'          ||       ||                               \n\
                   ''''     ''''    " + gappa_version() + ", (c) 2017-2018\n\
                                    by Lucas Czech and Pierre Barbera\n";
}

inline std::string gappa_title()
{
    return "gappa - Genesis Applications for Phylogenetic Placement Analysis";
}

inline std::string gappa_license()
{
    return "\
    gappa - Genesis Applications for Phylogenetic Placement Analysis.\n\
    Copyright (C) 2017-2018 Lucas Czech and HITS gGmbH\n\
    \n\
    This program is free software: you can redistribute it and/or modify\n\
    it under the terms of the GNU General Public License as published by\n\
    the Free Software Foundation, either version 3 of the License, or\n\
    (at your option) any later version.\n\
    \n\
    This program is distributed in the hope that it will be useful,\n\
    but WITHOUT ANY WARRANTY; without even the implied warranty of\n\
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the\n\
    GNU General Public License for more details.\n\
    \n\
    You should have received a copy of the GNU General Public License\n\
    along with this program.  If not, see <http://www.gnu.org/licenses/>.\n\
    \n\
    Contact:\n\
    Lucas Czech <lucas.czech@h-its.org>\n\
    Exelixis Lab, Heidelberg Institute for Theoretical Studies\n\
    Schloss-Wolfsbrunnenweg 35, D-69118 Heidelberg, Germany\n";
}

inline std::string gappa_citation()
{
    return "n/a";
}

#endif // include guard

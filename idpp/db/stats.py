"""
    idpp/db/stats.py

    Dylan Ross (dylan.ross@pnnl.gov)

    module with functions for computing database statistics to monitor 
    database growth and coverage, can also be run as a script
"""


import sys
from typing import Dict
import json

from matplotlib import pyplot as plt, patches as mpatches, rcParams
from matplotlib.axes import Axes as mplAxes

from idpp.db.util import IdPPdb


# set the font size and family for consistency with other plots
rcParams["font.size"] = 8
rcParams["font.family"] = "Roboto Condensed"


# define some colors
_C = {
    "none": "#000000",
    "rt": "#DD0000",
    "ccs": "#00DD00",
    "ms2": "#0000DD",
    "rt_ccs": "#DDDD00",
    "rt_ms2": "#DD00DD",
    "ccs_ms2": "#00DDDD",
    "rt_ccs_ms2": "#DDDDDD",
}


"""
    Bulk metrics
        - table sizes
        - other?
    potentially visualize them in some way 
    that reflects the database hierarchy 
"""


def compound_property_coverage(db: IdPPdb, include_none=True) -> Dict[str, int]:
    """
    plot the coverage of properties at the Compounds level

    ::

             |     X
             |     X                 X
         ^   |     X        X        X
         |   |     X     X  X  X     X
         #   |  X  X  X  X  X  X  X  X
             |  X  X  X  X  X  X  X  X
             |--------------------------
         RT  |  O  O     O  O
             |  |  |     |
        CCS  |  O  |  O  O     O
             |  |  |  |
        MS2  |  O  O  O           O   

                w  m  c  y  r  g  b  k

    Parameters
    ----------
    db : ``IdPPdb``
        IdPP database interface
    include_none : ``bool``
        include count of compounds without any associated properties
    
    Returns
    -------
    results : ``dict(str:int)``
        a dict mapping various property combinations to their
        corresponding counts
    """
    # various property combinations
    # build them from rules to reduce typing out a bunch of similar queries
    base = """--sqlite3
    SELECT
        COUNT(*)
    FROM
        (
            SELECT 
                cmpd_id,
                adduct_ids,
                adducts_,
                rt_ids,
                ccs_ids,
                ms2_ids
            FROM
                (
                    SELECT 
                        cmpd_id,
                        GROUP_CONCAT(adduct_id) AS adduct_ids,
                        GROUP_CONCAT(adduct) AS adducts_,
                        GROUP_CONCAT(rt_id) AS rt_ids,
                        GROUP_CONCAT(ccs_id) AS ccs_ids,
                        GROUP_CONCAT(ms2_id) AS ms2_ids
                    FROM 
                        Compounds 
                        JOIN 
                            Adducts USING(cmpd_id) 
                        LEFT JOIN 
                            RTs USING(adduct_id) 
                        LEFT JOIN 
                            CCSs USING(adduct_id) 
                        LEFT JOIN 
                            MS2Spectra USING(adduct_id) 
                    GROUP BY
                        cmpd_id
                )
            WHERE 
                rt_ids IS {not1} NULL
                AND ccs_ids IS {not2} NULL
                AND ms2_ids IS {not3} NULL
        )
    ;"""
    qrys = {
        "none": base.format(not1="", not2="", not3=""),
        "rt" : base.format(not1="NOT", not2="", not3=""),
        "ccs" : base.format(not1="", not2="NOT", not3=""),
        "ms2" : base.format(not1="", not2="", not3="NOT"),
        "rt_ccs" : base.format(not1="NOT", not2="NOT", not3=""),
        "rt_ms2" : base.format(not1="NOT", not2="", not3="NOT"),
        "ccs_ms2" : base.format(not1="", not2="NOT", not3="NOT"),
        "rt_ccs_ms2": base.format(not1="NOT", not2="NOT", not3="NOT"),
    }
    # run the queries
    results = {lbl: db.cur.execute(qry).fetchall()[0][0] for lbl, qry in qrys.items()}
    # make the plot
    fig, (ax_top, ax_bot) = plt.subplots(nrows=2, figsize=(2.5, 3.), height_ratios=(2, 1), sharex=True)
    # type annotations for convenience
    ax_top: mplAxes
    ax_bot: mplAxes
    w = 0.5
    r1, r2 = 0.6, 0.225
    #for i, (l, n) in enumerate(results.items()):
    groups = ["rt_ccs_ms2", "rt_ccs", "rt_ms2", "ccs_ms2", "rt", "ccs", "ms2", "none"]
    if not include_none:
        groups = groups[:-1]
    for i, l in enumerate(groups):
        n = results[l]
        ax_top.bar(i + 1, n, color=_C[l], edgecolor='k')
        if "rt" in l:
            ax_bot.add_artist(mpatches.Ellipse((i + 1, 0.8), 
                                               width=r1, height=r2, facecolor=_C[l], edgecolor="k"))
        if "ccs" in l:
            ax_bot.add_artist(mpatches.Ellipse((i + 1, 0.5), 
                                               width=r1, height=r2, facecolor=_C[l], edgecolor="k"))
        if "ms2" in l:
            ax_bot.add_artist(mpatches.Ellipse((i + 1, 0.2), 
                                               width=r1, height=r2, facecolor=_C[l], edgecolor="k"))
    for d in ['top', 'right']:
        ax_top.spines[d].set_visible(False)
    ax_top.set_xticks([_ + 1 for _ in range(8)])
    ax_top.set_xticklabels([])
    ax_top.ticklabel_format(axis='y', style='sci', scilimits=(0, 0))
    ax_top.set_ylabel('count')
    for d in ['top', 'right', 'bottom']:
        ax_bot.spines[d].set_visible(False)
    ax_bot.set_ylim([0, 1])
    ax_bot.set_yticks([0.2, 0.5, 0.8])
    ax_bot.set_yticklabels(['MS2', 'CCS', 'RT'])
    ax_bot.set_xticks([])
    plt.savefig(f'_figures/test_compound_property_coverage{"" if include_none else "_nonone"}.png', dpi=400, bbox_inches='tight')
    plt.tight_layout()
    return results
    

_get_colors = lambda lbl, n: [{"rt": lambda p: (p, 0.2, 0.2), "ccs": lambda p: (0.2, p, 0.2), "ms2": lambda p: (0.2, 0.2, p)}[lbl]((n - i) / (n)) for i in range(n)]
    

def property_source_distributions(db: IdPPdb) -> Dict[str, int]:
    """
    plot the number of entries in each property table over time
    based on timestamps in the tables

    ::

            +---------------------------+
        RTs |%%%|XXXXXXX|^v^v^v^v^v|0000|
            +---------------------------+   
                %  - Source 1
                X  - Source 2
                ^v - Source 3
                0  - Source 4

    Not exactly as in the example above, makes pie charts instead

    Parameters
    ----------
    db : ``IdPPdb``
        IdPP database interface

    Returns
    -------
    results : ``dict(...)``
        a dict mapping properties to lists of sources and counts of values
    """
    # set up a pattern for the query to select the sources counts
    qry = lambda prop: "SELECT src_name, COUNT(*) AS c FROM {tbl} JOIN Sources USING(src_id) GROUP BY src_name;".format(tbl={"rt":"RTs", "ccs":"CCSs", "ms2":"MS2Sources"}[prop])
    results = {
        "rt": db.cur.execute(qry("rt")).fetchall(),
        "ccs": db.cur.execute(qry("ccs")).fetchall(),
        "ms2": db.cur.execute(qry("ms2")).fetchall(),
    }
    # make the plots
    for prop in ["rt", "ccs", "ms2"]:
        fig, ax = plt.subplots(figsize=(2.5, 2.5))
        qdata = results[prop]
        if len(qdata) > 0:
            srcs, cnts = [], []
            for src, cnt in qdata:
                srcs.append(src)
                cnts.append(cnt)
            # sort by count
            srcs = [a for (a, b) in sorted(zip(srcs, cnts), key=lambda pair: pair[1], reverse=True)]
            cnts = sorted(cnts, reverse=True)
            n = len(cnts)
            # if there are more than 5 sources, sum together the rest into an "other" category
            if len(srcs) > 6:
                srcs = srcs[:5] + [f"other ({len(srcs[5:])})"]
                cnts = cnts[:5] + [sum(cnts[5:])]
                n = 6
            ax.pie(cnts, labels=srcs, colors=_get_colors(prop, n), explode=[0.1 + 0.05 * _ for _ in range(n)], wedgeprops={'linewidth': 1, 'ec': 'k'})
        ax.axis('off')
        plt.savefig(f"_figures/test_property_source_distributions_{prop}.png", dpi=400, bbox_inches='tight')
        plt.close()
    return results


# TODO: use fingerprinting and clustering or PCA or something to show the diversity of 
#       compounds in the database. This comparison can be made at the level of 
#       Compounds, Adducts, and property tables (RTs, CCSs, MS2Spectra)


def _main() -> None:

    # parse CLI args
    # <idpp_db> --gen_plots
    
    db = IdPPdb(sys.argv[1], read_only=True, enforce_idpp_ver=False)

    _ = compound_property_coverage(db, include_none=False)
    print(json.dumps(_, indent=2))

    _ = property_source_distributions(db)
    print(json.dumps(_, indent=2,))

    db.close()


if __name__ == "__main__":
    _main()

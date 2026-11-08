"""Plot suppression is local to an invocation and preserves numeric outputs."""
import click
from click.testing import CliRunner

from mlmm.core.plotting import add_plot_option, plot_enabled
from mlmm.io.plotly_image import write_plotly_image
from mlmm.io.trj2fig import save_outputs


def test_parent_policy_reaches_child_and_does_not_leak():
    states = []
    @click.command()
    @add_plot_option()
    def child():
        states.append(plot_enabled())
    @click.command()
    @add_plot_option()
    @click.pass_context
    def parent(ctx):
        states.append(plot_enabled())
        ctx.invoke(child)
    runner = CliRunner()
    assert runner.invoke(parent,["--no-plot"]).exit_code == 0
    assert states == [False,False]
    assert runner.invoke(parent).exit_code == 0
    assert states == [False,False,True,True]
    assert plot_enabled() is True


def test_no_plot_keeps_csv_without_rendering(tmp_path):
    class ForbiddenFigure:
        def write_html(self,*args,**kwargs):
            raise AssertionError("HTML rendering should be skipped")
        def to_json(self,*args,**kwargs):
            raise AssertionError("Chrome export should be skipped before serialization")
    csv = tmp_path/"energies.csv"
    html,png = tmp_path/"plot.html",tmp_path/"plot.png"
    @click.command()
    @add_plot_option()
    def command():
        assert write_plotly_image(ForbiddenFigure(),png) is None
        save_outputs([csv,html,png],ForbiddenFigure(),[1.,2.],[0.,1.],"hartree",True)
    result = CliRunner().invoke(command,["--no-plot"])
    assert result.exit_code == 0, result.exception
    assert len(csv.read_text().splitlines()) == 3
    assert not html.exists() and not png.exists()


def test_all_dispatch_to_a_fresh_click_context_inherits_no_plot():
    from mlmm.workflows.all import _run_cli_main
    states = []
    @click.command()
    @add_plot_option()
    def child():
        states.append(plot_enabled())
    @click.command()
    @add_plot_option()
    def parent():
        assert _run_cli_main("child",child,[]) == 0
    result = CliRunner().invoke(parent,["--no-plot"])
    assert result.exit_code == 0, result.exception
    assert states == [False]


def test_energy_level_without_diagram_keeps_numerical_outputs():
    from mlmm.workflows._all_helpers import build_energy_level_dict

    result = build_energy_level_dict(labels=["R", "TS", "P"],
        energies_au=[-1.0, -0.9, -1.1], ref_energy=-1.0, au_to_kcal=627.5,
        diagram_path=None, structures={"R": "r.pdb", "TS": "ts.pdb", "P": "p.pdb"})
    assert result["diagram"] is None
    assert result["energies_au"] == [-1.0, -0.9, -1.1]
    assert result["barrier_kcal"] > 0


def test_root_promotion_keeps_omitted_image_metadata(tmp_path):
    from mlmm.workflows._all_helpers import promote_diag_for_root

    original = {"name": "energy_diagram_MLIP", "image": None, "image_written": False,
                "energies_au": [-1.0, -0.9, -1.1]}
    promoted = promote_diag_for_root(original, "energy_diagram_MLIP", tmp_path)
    assert promoted["image"] is None
    assert promoted["image_written"] is False
    assert promoted["energies_au"] == original["energies_au"]
    assert original["name"] == "energy_diagram_MLIP"

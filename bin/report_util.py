
from dominate.tags import div, h1, h4, p, code, hr
from dominate.tags import figure, img, a
import base64
from ezcharts.layout.snippets import DataTable, Grid, Tabs
from ezcharts.layout.resource import (ImageResource, ScriptResource, StyleResource, ThemeResource, Resource)
from ezcharts.layout.util import inline, load_json, transpile, cls, css
from ezcharts.layout.base import IClasses, IStyles, Snippet
from ezcharts.layout.snippets.banner import IBannerClasses, IBannerStyles, IBadgeClasses, IBadgeStyles, Badge
from ezcharts.components.reports.labs import ILabsAddendumClasses
from pathlib import Path
import json
from PIL import Image



class EIBanner(Snippet):
    """A styled div tag containing a heading and badges."""

    TAG = 'div'

    def __init__(
        self,
        report_title: str,
        workflow_name: str,
        styles: IBannerStyles = IBannerStyles(),
        classes: IBannerClasses = IBannerClasses(),
        default_content: bool = True
    ) -> None:
        """Create styled banner."""
        super().__init__(
            styles=styles,
            classes=classes,
            className=classes.container,
            style=styles.container)

        with self:
            with div(className=classes.inner, style=styles.inner):
                if not default_content:
                    return
                h1(report_title)
                p(
                    f"Results generated through the {workflow_name.upper()} nextflow "
                    "workflow developed by Earlham Institute.",
                    className="py-3 fs-5")
                self.badges = div(className="d-flex flex-wrap")

    def add_badge(
        self,
        title: str,
        bg_class=None,
        styles: IBadgeStyles = IBadgeStyles(),
        classes: IBadgeClasses = IBadgeClasses(),
    ) -> None:
        """Add a badge to the banner."""
        with self.badges:
            Badge(title, styles=styles, classes=classes, bg_class=bg_class)


class EILabsAddendum(Snippet):
    """A styled footer component for use in a Report."""

    TAG = 'div'

    def __init__(
        self,
        workflow_name: str,
        workflow_version: str,
        classes: ILabsAddendumClasses = ILabsAddendumClasses(),
        use_defaults: bool = True
    ) -> None:
        """Create tag."""
        super().__init__(
            styles=None,
            classes=classes,
            className=classes.container)

        with self:
            self.container = div(className=classes.inner)

        if use_defaults:
            with self.container:
                h4('About this report', className="pb-3")
                p(
                    "This report was produced using the ",
                    code(f"EarlhamInst/{workflow_name}"),
                    f" nextflow workflow (v{workflow_version})."
                )
                p(
                    "The analysis restuls are not "
                    "intended for use for health assessment or to "
                    "diagnose, treat, mitigate, cure or prevent any "
                    "disease or condition.")



def plots_from_image_files(path, meta=None, ncol=1, suffix=['*.png'], widths=None):
    """Create a plots section which uses image files.
    PARAMS:
        meta: can be sample or group, by which images are shown in tabs
        ncol: the number of columns for Grid
        suffix: list of suffix used to match the image files
        widths: the column withs of Grid columns, e.g.['500','700']
    """
    if not widths:
        widths = [1200//ncol]*ncol

    # if meta in ['sample', 'group']:
    if meta:
        tabs = Tabs()
        for folder in sorted([f for f in path.iterdir() if f.is_dir()]):
            if folder.name.startswith(f'{meta}_'):
                sample_id = folder.name.split('_', 1)[1]
                with tabs.add_tab(sample_id):
                    # pathes = [next(folder.glob('*'+sf)) for sf in suffix]
                    pathes = sorted([file for sf in suffix for file in folder.glob(sf) if file.is_file()])
                    Np = len(pathes)
                    widths = (widths*(Np//len(widths)+1))[0:Np]
                    for i in [r*ncol for r in range(Np//ncol+Np%ncol)]:
                        with Grid(columns=ncol if Np//ncol>0 else Np):
                            pathes_r = [pathes[i+s] for s in range(ncol) if i+s < Np]
                            widths_r = [widths[i+s] for s in range(ncol) if i+s < Np]
                            for img_path, width in zip(pathes_r, widths_r):
                                with Image.open(img_path) as image:
                                    width_p, height_p = image.size
                                    width = min(int(width), width_p)
                                with open(img_path, 'rb') as fh:
                                    b64img = base64.b64encode(fh.read()).decode()
                                    with figure(cls='text-center'):
                                        img(src=f'data:image/png;base64,{b64img}', width=width)
                                        # a(img(src=f'data:image/png;base64,{b64img}', width=width), href=f'{img_path}')
    else:
        pathes = sorted([file for sf in suffix for file in path.glob(sf) if file.is_file()])
        Np = len(pathes)
        widths = (widths*(Np//len(widths)+1))[0:Np]
        for i in [r*ncol for r in range(Np//ncol+Np%ncol)]:        
            with Grid(columns=ncol if Np//ncol>0 else Np):
                pathes_r = [pathes[i+s] for s in range(ncol) if i+s < Np]
                widths_r = [widths[i+s] for s in range(ncol) if i+s < Np]
                for img_path, width in zip(pathes_r, widths_r):
                    with Image.open(img_path) as image:
                        width_p, height_p = image.size
                        width = min(int(width), width_p)                    
                    with open(img_path, 'rb') as fh:
                        b64img = base64.b64encode(fh.read()).decode()
                        with figure(cls='text-center'):
                            img(src=f'data:image/png;base64,{b64img}', width=width)
                            # a(img(src=f'data:image/png;base64,{b64img}', width=width), href=f'{img_path}')        


class EILogo(div):
    """EI logo element."""

    def __init__(self) -> None:
        """Create a div with an SVG logo inside."""
        super().__init__(
            # inline(ImageResource('EI_logo.svg').data_file),
            # inline(Resource(path='images/EI_logo.svg', loader=inline)),
            inline(Path('images/EI_logo.svg')),
            # inline(Path('/home/wuhh/nf-core-eisca/bin/images/EI_logo.svg')),
            tagname='div',
            style="width: 35px; height: 35px;",
            className="d-flex",
            alt="EI Logo")
        


def show_analysis_parameters(params_file):
    with open(params_file, 'r') as file:
        # params = load_json(file)
        params = json.load(file)
        params = ' | '.join([f"{k} {v}" for k, v in params.items()])
        # h4('About this report', className="pb-3")
        p(
            hr(),
            "The analysis parameters: ", 
            code(params)
        )

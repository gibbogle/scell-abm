#include "field.h"
#include "log.h"

#include "global.h"

LOG_USE();

//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------
//Field::Field(QWidget *page2D)
Field::Field(QWidget *aParent) : QWidget(aParent)
{
//	field_page = page2D;
    field_page = aParent;
    axis = Z_AXIS;
    fraction = 0;
    constituent = OXYGEN;
    slice_changed = true;
    setConcPlot(false);
    setVolPlot(false);
    setOxyPlot(false);
    pGconc = NULL;
    pGvol = NULL;
    pGoxy = NULL;
    ifield = 0;
    view = new MyQGraphicsView(field_page);
    constituent_rb_list = NULL;
    vbox_constituent = NULL;
    buttonGroup_constituent = new QButtonGroup;
    data = NULL;
}

//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------
Field::~Field()
{
}

//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------
bool Field::isConcPlot()
{
    return useConcPlot;
}

//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------
void Field::setConcPlot(bool status)
{
    useConcPlot = status;
}

//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------
bool Field::isVolPlot()
{
    return useVolPlot;
}

//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------
void Field::setVolPlot(bool status)
{
    useVolPlot = status;
}

//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------
bool Field::isOxyPlot()
{
    return useOxyPlot;
}

//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------
void Field::setOxyPlot(bool status)
{
    useOxyPlot = status;
}

//------------------------------------------------------------------------------------------------
// To create the group of radiobuttons for constituent selection.
// This uses information about active constituents fetched from the DLL.
//------------------------------------------------------------------------------------------------
void Field::setConstituentButtons(QGroupBox *gbox, QButtonGroup *bg, QVBoxLayout **vbox, QRadioButton ***rb_list, QString tag)
{
    int ivar;
    QString name, str;
    int **ip;
    QRadioButton **p;
    QRadioButton *rb;

    p = *rb_list;
    LOG_QMSG("setConstituentButtons: " + tag);
    if (p) {
        LOG_MSG("rb_list not NULL, delete it");
        for (ivar=0; ivar<Global::nvars_used; ivar++) {
            rb = p[ivar];
            bg->removeButton(rb);
            delete rb;
        }
        delete p;
    }
    if (!*vbox) {
        LOG_MSG("vbox = NULL, create it");
        *vbox = new QVBoxLayout;
        gbox->setLayout(*vbox);
    }
    name = "rb_constituent_"+tag;
    LOG_QMSG(name);
    *rb_list = new QRadioButton*[Global::nvars_used];
    p = *rb_list;
//    sprintf(msg,"rb_list: %p vbox: %p bg: %p nvars_used: %d",p,*vbox,bg,Global::nvars_used);
//    LOG_MSG(msg);
    for (ivar=0; ivar<Global::nvars_used; ivar++) {
        str = Global::var_string[ivar];
        p[ivar] = new QRadioButton;
        p[ivar]->setText(str);
        p[ivar]->setObjectName(name+ivar);
        (*vbox)->addWidget(p[ivar]);
        p[ivar]->setEnabled(true);
        bg->addButton(p[ivar],ivar);
    }
    p[1]->setChecked(true);   // Oxygen
    QRect rect = gbox->geometry();
    rect.setHeight(25*Global::nvars_used);
    gbox->setGeometry(rect);
    gbox->show();
}

//------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------
void Field::selectConstituent()
{
    int iconst, res;
    QStringList items;

    LOG_MSG("selectConstituent");
    for (iconst=0; iconst<Global::nvars_used; iconst++) {
        if (iconst == constituent) continue;
        items << Global::var_string[iconst];
    }
    bool ok;
    QString item = QInputDialog::getItem(this, tr("QInputDialog::getItem()"),
                                          tr("Constituent:"), items, 0, false, &ok);
    if (ok && !item.isEmpty()) {
        for (iconst=0; iconst<Global::nvars_used; iconst++) {
            if (item == Global::var_string[iconst]) {
                constituent = iconst;
//                if (useConcPlot)
//                    updateConcPlot();
//                break;
            }
        }
    }
}


//------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------
void Field::setConstituent(QAbstractButton *button)
{
    int res;
//    QMessageBox msgBox;
//    msgBox.setText("setConstituent");
//    msgBox.exec();
    LOG_MSG("setConstituent");
    int prev_constituent = constituent;
    QString text = button->text();
    for (int ivar=0; ivar<Global::nvars_used; ivar++) {
        if (text == Global::var_string[ivar]) {
            constituent = ivar;
        }
    }

    if (constituent != prev_constituent) {
        constituent_changed = true;
        LOG_MSG("setConstituent");
        displayField(hour,&res);
    }
}

//------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------
void Field::setPlane(QAbstractButton *button)
{
    int res;
//    QMessageBox msgBox;
//    msgBox.setText("setPlane");
//    msgBox.exec();
    LOG_MSG("setPlane");
    QString text = button->text();
	int prev_axis = axis;
    if (text.compare("X-Y") == 0)
        axis = Z_AXIS;
    else if (text.compare("Y-Z") == 0)
        axis = X_AXIS;
    else if (text.compare("X-Z") == 0)
        axis = Y_AXIS;
	if (axis != prev_axis) {
        slice_changed = true;
        LOG_MSG("setPlane");
        displayField(hour,&res);
	}
}

//------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------
void Field::setFraction(QString text)
{
    int res;
	double prev_fraction = fraction;
	fraction = text.toDouble();
	if (fraction != prev_fraction) {
        displayField(hour,&res);
	}
}

//------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------
void Field::setSliceChanged()
{
    slice_changed = true;
}

//------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------
void Field::setSaveImages(bool save)
{
    save_images = save;
}

//------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------
void Field::setUseLogScale(bool use_logscale)
{
    use_log = use_logscale;
}

//-----------------------------------------------------------------------------------------
// New version, site/cell size is fixed, the blob grows
//-----------------------------------------------------------------------------------------
void Field::displayField(int hr, int *res)
{
    QGraphicsScene* scene = new QGraphicsScene(QRect(0, 0, CANVAS_WIDTH, CANVAS_WIDTH));
    QBrush brush;
    int i, xindex, yindex, ix, iy, w, rgbcol[3], ichemo;
    double xp, yp, d0, d, volume, scale, cmin, cmax, rmax;
    double a, b, Wc;
    int Nc;

    ichemo = Global::GUI_to_DLL_index[constituent];
    LOG_QMSG("displayField: " + QString::number(constituent) + "-->" + QString::number(ichemo));
    use_log = false;    // temporary
    *res = 0;
    hour = hr;
	if (slice_changed) {
//        get_fieldinfo(&NX, &axis, &fraction, &nsites, &nconst, const_used, res);    // Note: const_used[] is no longer used
        if (*res != 0) {
            printf("Error: get_fieldinfo: FAILED\n");
            LOG_MSG("Error: get_fieldinfo: FAILED");
            return;
        }
        if (nconst != MAX_CONC) {
            sprintf(msg,"Error: displayField: get_fieldinfo: MAX_CONC != MAX_CHEMO: %d %d",MAX_CONC,Global::MAX_CHEMO);
            LOG_MSG(msg);
            exit(1);
        }
        if (NEXTRA != Global::N_EXTRA) {
            sprintf(msg,"Error: displayField: NEXTRA != N_EXTRA: %d %d",NEXTRA,Global::N_EXTRA);
            LOG_MSG(msg);
            exit(1);
        }
        if (this->data) {
            free(this->data);
        }
        this->data = (FIELD_DATA *)malloc(nsites*sizeof(FIELD_DATA));
//        get_fielddata(&axis, &fraction, &nsites, &nconst, this->data, res);
        if (*res != 0) {
            LOG_MSG("Error: get_fielddata: FAILED");
            return;
        }
        slice_changed = false;
    }

    if (axis == X_AXIS) {           // Y-Z plane
        xindex = 1;
        yindex = 2;
    } else if (axis == Y_AXIS) {   // X-Z plane
        xindex = 0;
        yindex = 2;
    } else if (axis == Z_AXIS) {   // X-Y plane
        xindex = 0;
        yindex = 1;
    }

/*
 NX = size of lattice
 Nc = # of sites to fill the canvas from side to side (or top to bottom) = (2/3)NX
 Wc = canvas width (pixels)
 w = site width = Wc/Nc
 xp = a.ix + b
 yp = a.iy + b
 blob centre at (NX/2,NX/2) maps to canvas centre at (Wc/2,Wc/2)
 => Wc/2 = a.NX/2 + b
 The width of Nc sites maps to the canvas width
 => Wc = a.Nc
 => a = Wc/Nc, b = Wc/2 - a.NX/2
*/
    Nc = (2*NX)/3;
    Wc = CANVAS_WIDTH;
    w = Wc/Nc;
    a = w;
    b = Wc/2 - a*NX/2;
    d0 = w*Global::dfraction;
    cmin = 1.0e10;
    cmax = 0;
    rmax = 0;
    for (i=0; i<nsites; i++) {
        rmax = MAX(rmax,data[i].conc[GROWTH_RATE]);     // GROWTH_RATE = Global::MAX_CHEMO+1
        cmin = MIN(MAX(cmin,0),data[i].conc[ichemo]);
        cmax = MAX(cmax,data[i].conc[ichemo]);
        // Flip it
//        data[i].site[yindex] = NX - data[i].site[yindex];
    }
    brush.setStyle(Qt::SolidPattern);
    brush.setColor(QColor(0,0,0));
    scene->addRect(0,0,CANVAS_WIDTH,CANVAS_WIDTH,Qt::NoPen, brush);
    view->setScene(scene);
    view->setGeometry(QRect(0, 0, 700, 700));
    if (cmax == 0) {
        view->show();
        return;
    }
    for (i=0; i<nsites; i++) {
        ix = this->data[i].site[xindex];
        iy = this->data[i].site[yindex];
        xp = int(a*ix + b - w);
        yp = int(a*iy + b - w);
        chooseFieldColor(data[i].conc[ichemo],cmin,cmax,use_log,rgbcol);
//        sprintf(msg,"c: %f %f %f rgbcol: %d %d %d",data[i].conc[constituent],cmin,cmax,rgbcol[0],rgbcol[1],rgbcol[2]);
//        LOG_MSG(msg);
        brush.setColor(QColor(rgbcol[0],rgbcol[1],rgbcol[2]));
        scene->addRect(xp,yp,w,w,Qt::NoPen, brush);
    }
    for (i=0; i<nsites; i++) {
        ix = this->data[i].site[xindex];
        iy = this->data[i].site[yindex];
        xp = int(a*ix + b - w);
        yp = int(a*iy + b - w);
        volume = this->data[i].volume;      // = 0 if there is no cell
        if (volume > 0) {
            scale = pow(volume,0.3333);
            d = scale*d0;   // fix this - need to change d0
            double f = data[i].conc[GROWTH_RATE]/rmax;
            chooseRateColor(f,rgbcol);
            brush.setColor(QColor(rgbcol[0],rgbcol[1],rgbcol[2]));
            scene->addEllipse(xp+(w-d)/2,yp+(w-d)/2,d,d,Qt::NoPen, brush);
        }
    }
    view->show();
    if (save_images) {
        scene->clearSelection();                                                  // Selections would also render to the file
        scene->setSceneRect(scene->itemsBoundingRect());                          // Re-shrink the scene to it's bounding contents
        QImage image(scene->sceneRect().size().toSize(), QImage::Format_ARGB32);  // Create the image with the exact size of the shrunk scene
        image.fill(Qt::transparent);                                              // Start all pixels transparent

        QPainter painter(&image);
        scene->render(&painter);
        ifield++;
        char filename[] = "image/field0000.png";
        char numstr[5];
        sprintf(numstr,"%04d",hour);
        for (int i=0; i<4; i++)
            filename[11+i] = numstr[i];
        image.save(filename);
    }
}


//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------
void Field::chooseFieldColor(double c, double cmin, double cmax, bool use_logscale, int rgbcol[])
{
    double f, denom, logcmin, logc;
    int rgb_lo[3], rgb_hi[3], i;

    if (use_logscale) {
        if (cmin == cmax) {
            f = 1;
        } else {
//            if (cmin > 0.0001)
//                logcmin = log(cmin);
//            else
//                logcmin = 1.0e6;
//            if (c > 0.0001)
//                logc = log(c);
//            else
//                logc = 1.0e6;
            cmin = max(cmin, 0.00001);
            c = max(c, 0.00001);
            denom = (log(cmax) - log(cmin));
            if (denom < 0.001)
                f = 1;
            else
                f = (log(c) - log(cmin))/denom;
        }
    } else {
        f = c/cmax;
    }
    if (constituent == OXYGEN) {
        rgb_hi[0] =   0; rgb_hi[1] =   0; rgb_hi[2] = 0;
        rgb_lo[0] = 255; rgb_lo[1] =   0; rgb_lo[2] = 0;
        for (i=0; i<3; i++) {
            rgbcol[i] = int((1-f)*rgb_lo[i] + f*rgb_hi[i]);
            if (rgbcol[i] < 0 || rgbcol[i] > 255) {
                sprintf(msg,"chooseFieldColor: %f %f %f %f %d %d",c,cmin,cmax,f,i,rgbcol[i]);
                LOG_MSG(msg);
                exit(1);
            }
        }
    } else {
        rgb_hi[0] =   0; rgb_hi[1] =   255; rgb_hi[2] = 255;
        rgb_lo[0] =   0; rgb_lo[1] =   0; rgb_lo[2] = 0;
        for (i=0; i<3; i++) {
            rgbcol[i] = int((1-f)*rgb_lo[i] + f*rgb_hi[i]);
            if (rgbcol[i] < 0 || rgbcol[i] > 255) {
                sprintf(msg,"chooseFieldColor: %f %f %f %f %d %d",c,cmin,cmax,f,i,rgbcol[i]);
                LOG_MSG(msg);
                exit(1);
            }
        }
    }
}

//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------
void Field::chooseRateColor(double f, int rgbcol[])
{
    int rgb_lo[3], rgb_hi[3], i;

    rgb_hi[0] = 0; rgb_hi[1] = 255; rgb_hi[2] = 0;
    rgb_lo[0] = 0; rgb_lo[1] =  64; rgb_lo[2] = 0;
    for (i=0; i<3; i++) {
        rgbcol[i] = int((1-f)*rgb_lo[i] + f*rgb_hi[i]);
    }
}

//-----------------------------------------------------------------------------------------
// Now 'constituent' is an index of the active constituents: 0 - nvars_used
//-----------------------------------------------------------------------------------------
void Field::getTitle(int iconst, QString *title)
{
    QString name = Global::var_string[iconst];
    if (Global::GUI_to_DLL_index[iconst] <= Global::MAX_CHEMO) {
        *title = name + " Concentration";
    } else {
        *title = name;
    }
}

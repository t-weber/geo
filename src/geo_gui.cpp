/**
 * convex hull test program
 * @author Tobias Weber (orcid: 0000-0002-7230-1932)
 * @date 15-Aug-2020
 * @license: see 'LICENSE' file
 */

#include "geo_gui.h"

#include <QApplication>
#include <QMenuBar>
#include <QLabel>
#include <QStatusBar>
#include <QMouseEvent>
#include <QFileDialog>
#include <QMessageBox>
#include <QSettings>

#include <locale>
#include <memory>
#include <array>
#include <vector>
#include <iostream>

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>

#include "geo_algos.h"


using t_real = double;
using t_vec = m::vec<t_real, std::vector>;
using t_mat = m::mat<t_real, std::vector>;

namespace ptree = boost::property_tree;

const t_real g_eps = 1e-5;


// ----------------------------------------------------------------------------
// #define HULL_CHECK

#ifdef HULL_CHECK
static double side_of_line(const QLineF& line, const QPointF& pt)
{
	QPointF dir1 = line.p2() - line.p1();
	QPointF dir2 = pt - line.p1();

	return dir1.x()*dir2.y() - dir1.y()*dir2.x();
}


static bool all_points_on_same_side(const QLineF& line, const std::vector<QPointF>& hullvertices)
{
	// find a reference vertex which is sufficiently far from the line
	std::optional<double> side;
	for(const QPointF& vert : hullvertices)
	{
		if(!side)
		{
			double curside = side_of_line(line, vert);
			if(std::abs(curside) > g_eps)
				side = curside;
		}

		if(side)
			break;
	}

	if(!side)
		return true;


	// are all other vertices on the same side as the reference vertex (or on the line)?
	for(const QPointF& vert : hullvertices)
	{
		double curside = side_of_line(line, vert);
		if(std::signbit(*side) != std::signbit(curside) && std::abs(curside) > g_eps)
			return false;
	}

	return true;
}
#endif
// ----------------------------------------------------------------------------



// ----------------------------------------------------------------------------

Vertex::Vertex(const QPointF& pos, double rad) : m_rad{rad}
{
	setPos(pos);
	setFlags(flags() | QGraphicsItem::ItemIsMovable | QGraphicsItem::ItemIsSelectable);
}


Vertex::~Vertex()
{
}


QRectF Vertex::boundingRect() const
{
	return QRectF{-m_rad/2., -m_rad/2., m_rad, m_rad};
}


void Vertex::paint(QPainter* painter, const QStyleOptionGraphicsItem*, QWidget*)
{
	std::array<QColor, 2> colours =
	{
		QColor::fromRgbF(0.,0.,1.),
		QColor::fromRgbF(0.,0.,0.),
	};

	QRadialGradient grad{};
	grad.setCenter(0., 0.);
	grad.setRadius(m_rad);

	for(std::size_t col=0; col<colours.size(); ++col)
		grad.setColorAt(col/double(colours.size()-1), colours[col]);

	painter->setBrush(grad);
	painter->setPen(*colours.rbegin());

	painter->drawEllipse(-m_rad/2., -m_rad/2., m_rad, m_rad);
}

// ----------------------------------------------------------------------------




// ----------------------------------------------------------------------------

HullView::HullView(QGraphicsScene *scene, QWidget *parent) : QGraphicsView(scene, parent),
	m_scene{scene}
{
	setHorizontalScrollBarPolicy(Qt::ScrollBarAlwaysOn);
	setVerticalScrollBarPolicy(Qt::ScrollBarAlwaysOn);

	setInteractive(true);
	setMouseTracking(true);

	setBackgroundBrush(QBrush{QColor::fromRgbF(0.95, 0.95, 0.95, 1.)});
}


HullView::~HullView()
{
}


void HullView::resizeEvent(QResizeEvent *evt)
{
	QPointF pt1{mapToScene(QPoint{0,0})};
	QPointF pt2{mapToScene(QPoint{evt->size().width(), evt->size().height()})};

	const double padding = 16;

	// include bounds given by vertices
	for(const Vertex* vertex : m_vertices)
	{
		QPointF vertexpos = vertex->scenePos();

		if(vertexpos.x() < pt1.x())
			pt1.setX(vertexpos.x() -  padding);
		if(vertexpos.x() > pt2.x())
			pt2.setX(vertexpos.x() +  padding);
		if(vertexpos.y() < pt1.y())
			pt1.setY(vertexpos.y() -  padding);
		if(vertexpos.y() > pt2.y())
			pt2.setY(vertexpos.y() +  padding);
	}

	setSceneRect(QRectF{pt1, pt2});
}



void HullView::AddVertex(const QPointF& pos)
{
	Vertex *vertex = new Vertex{pos};
	m_vertices.insert(vertex);
	m_scene->addItem(vertex);
}


void HullView::mousePressEvent(QMouseEvent *evt)
{
	QPoint posVP = evt->pos();
	QPointF posScene = mapToScene(posVP);

	QList<QGraphicsItem*> items = this->items(posVP);
	QGraphicsItem* item = nullptr;
	bool item_is_vertex = false;

	for(int itemidx=0; itemidx<items.size(); ++itemidx)
	{
		item = items[itemidx];
		item_is_vertex = m_vertices.find(static_cast<Vertex*>(item)) != m_vertices.end();
		if(item_is_vertex)
			break;
	}

	// only select vertices
	if(!item_is_vertex)
		item = nullptr;


	if(evt->button() == Qt::LeftButton)
	{
		// if no vertex is at this position, create a new one
		if(!item)
		{
			AddVertex(posScene);
			m_dragging = true;
			UpdateAll();
		}

		else
		{
			// vertex is being dragged
			if(item_is_vertex)
			{
				m_dragging = true;
			}
		}
	}
	else if(evt->button() == Qt::RightButton)
	{
		// if a vertex is at this position, remove it
		if(item && item_is_vertex)
		{
			m_scene->removeItem(item);
			m_vertices.erase(static_cast<Vertex*>(item));
			delete item;
			UpdateAll();
		}
	}

	QGraphicsView::mousePressEvent(evt);
}


void HullView::mouseReleaseEvent(QMouseEvent *evt)
{
	if(evt->button() == Qt::LeftButton)
		m_dragging = false;

	UpdateAll();
	QGraphicsView::mouseReleaseEvent(evt);
}


void HullView::mouseMoveEvent(QMouseEvent *evt)
{
	QGraphicsView::mouseMoveEvent(evt);

	if(m_dragging)
	{
		QResizeEvent evt{size(), size()};
		resizeEvent(&evt);
		UpdateAll();
	}

	QPoint posVP = evt->pos();
	QPointF posScene = mapToScene(posVP);
	emit SignalMouseCoordinates(posScene.x(), posScene.y());
}


void HullView::SetCalculateHull(bool b)
{
	m_calchull = b;
	UpdateHull();
}


void HullView::SetCalculateVoronoiVertices(bool b)
{
	m_calcvoronoivertices = b;
	UpdateDelaunay();
}


void HullView::SetCalculateVoronoiRegions(bool b)
{
	m_calcvoronoiregions = b;
	UpdateDelaunay();
}


void HullView::SetCalculateDelaunay(bool b)
{
	m_calcdelaunay = b;
	UpdateDelaunay();
}


void HullView::SetHullCalculationMethod(HullCalculationMethod m)
{
	m_hullcalculationmethod = m;
	UpdateHull();
}


void HullView::SetDelaunayCalculationMethod(DelaunayCalculationMethod m)
{
	m_delaunaycalculationmethod = m;
	UpdateDelaunay();
}


void HullView::ClearVertices()
{
	for(Vertex* vertex : m_vertices)
	{
		m_scene->removeItem(vertex);
		delete vertex;
	}
	m_vertices.clear();

	UpdateAll();
}


void HullView::UpdateAll()
{
	UpdateDelaunay();
	UpdateHull();
}


void HullView::UpdateHull()
{
	// remove previous hull
	for(QGraphicsItem* hullItem : m_hull)
	{
		m_scene->removeItem(hullItem);
		delete hullItem;
	}
	m_hull.clear();

	if(!m_calchull || m_vertices.size() < 3)
		return;

	std::vector<t_vec> vertices;
	vertices.reserve(m_vertices.size());
	std::transform(m_vertices.begin(), m_vertices.end(), std::back_inserter(vertices),
		[](const Vertex* vert) -> t_vec { return m::create<t_vec>({vert->x(), vert->y()}); } );


	std::vector<std::vector<t_vec>> hull;

	switch(m_hullcalculationmethod)
	{
		case HullCalculationMethod::QHULL:
			std::tie(std::ignore, hull, std::ignore) = calc_delaunay<t_vec>(2, vertices, true);
			break;
		case HullCalculationMethod::CONTOUR:
			hull.emplace_back(calc_hull_contour<t_vec>(vertices));
			break;
		case HullCalculationMethod::ITERATIVE:
			//hull.emplace_back(calc_hull_iterative<t_vec>(vertices));
			hull.emplace_back(calc_hull_iterative_bintree<t_vec>(vertices));
			break;
		case HullCalculationMethod::RECURSIVE:
			hull.emplace_back(calc_hull_recursive<t_vec>(vertices));
			break;
		default:
			QMessageBox::critical(this, "Error", "Unknown hull calculation method.");
			break;
	}


#ifdef HULL_CHECK
	std::vector<QPointF> hullvertices;
	for(const auto& thetriag : hull)
		for(std::size_t idx1=0; idx1<thetriag.size(); ++idx1)
			hullvertices.emplace_back(QPointF{thetriag[idx1][0], thetriag[idx1][1]});
#endif

	// convex hull
	QPen penHull;
	penHull.setWidthF(2.);

	for(const auto& thetriag : hull)
	{
		for(std::size_t idx1=0; idx1<thetriag.size(); ++idx1)
		{
			std::size_t idx2 = idx1+1;
			if(idx2 >= thetriag.size())
				idx2 = 0;
			if(idx1 == idx2)
				continue;

			QLineF line{QPointF{thetriag[idx1][0], thetriag[idx1][1]}, QPointF{thetriag[idx2][0], thetriag[idx2][1]}};
#ifdef HULL_CHECK
			if(!all_points_on_same_side(line, hullvertices))
				continue;
#endif

			QGraphicsItem *item = m_scene->addLine(line, penHull);
			m_hull.insert(item);
		}
	}
}


void HullView::UpdateDelaunay()
{
	// remove previous triangulation
	for(QGraphicsItem* item : m_delaunay)
	{
		m_scene->removeItem(item);
		delete item;
	}
	m_delaunay.clear();

	// remove previous voronoi vertices
	for(QGraphicsItem* item : m_voronoi)
	{
		m_scene->removeItem(item);
		delete item;
	}
	m_voronoi.clear();


	if((!m_calcdelaunay && !m_calcvoronoivertices && !m_calcvoronoiregions) || m_vertices.size() < 4)
		return;

	std::vector<t_vec> vertices;
	vertices.reserve(m_vertices.size());
	std::transform(m_vertices.begin(), m_vertices.end(), std::back_inserter(vertices),
		[](const Vertex* vert) -> t_vec { return m::create<t_vec>({vert->x(), vert->y()}); } );


	std::vector<t_vec> voronoi{};
	std::vector<std::vector<t_vec>> triags{};
	std::vector<std::set<std::size_t>> neighbours{};

	switch(m_delaunaycalculationmethod)
	{
		case DelaunayCalculationMethod::QHULL:
			std::tie(voronoi, triags, neighbours) = calc_delaunay<t_vec>(2, vertices, false);
			break;
		case DelaunayCalculationMethod::PARABOLIC:
			std::tie(voronoi, triags, neighbours) = calc_delaunay_parabolic<t_vec>(vertices);
			break;
		default:
			QMessageBox::critical(this, "Error", "Unknown Delaunay calculation method.");
			break;
	}


	const t_real itemRad = 7.;

	if(m_calcvoronoivertices)
	{
		QPen penVoronoi;
		penVoronoi.setStyle(Qt::SolidLine);
		penVoronoi.setWidthF(1.);

		QPen penCircle;
		penCircle.setStyle(Qt::DotLine);
		penCircle.setWidthF(1.);
		penCircle.setColor(QColor::fromRgbF(1.,0.,0.));

		QBrush brushVoronoi;
		brushVoronoi.setStyle(Qt::SolidPattern);
		brushVoronoi.setColor(QColor::fromRgbF(1.,0.,0.));

		// voronoi vertices
		for(std::size_t idx=0; idx<voronoi.size(); ++idx)
		{
			const t_vec& voronoivert = voronoi[idx];

			QPointF voronoipt{voronoivert[0], voronoivert[1]};
			QGraphicsItem *voronoiItem = m_scene->addEllipse(
				voronoipt.x()-itemRad/2., voronoipt.y()-itemRad/2., itemRad, itemRad, penVoronoi, brushVoronoi);
			m_voronoi.insert(voronoiItem);

			// circles
			if(idx < triags.size())
			{
				const auto& triag = triags[idx];
				if(triag.size() >= 3)
				{
					t_real rad = m::norm(voronoivert-triag[0]);

					QGraphicsItem *voronoiCircle = m_scene->addEllipse(
						voronoipt.x()-rad, voronoipt.y()-rad, rad*2., rad*2., penCircle);
					m_voronoi.insert(voronoiCircle);
				}
			}
		}
	}


	if(m_calcvoronoiregions && neighbours.size()==voronoi.size())
	{
		QPen penVoronoi;
		penVoronoi.setStyle(Qt::SolidLine);
		penVoronoi.setWidthF(1.);
		penVoronoi.setColor(QColor::fromRgbF(1.,0.,0.));

		QPen penVoronoiUnbound;
		penVoronoiUnbound.setStyle(Qt::DashLine);
		penVoronoiUnbound.setWidthF(1.);
		penVoronoiUnbound.setColor(QColor::fromRgbF(1.,0.,0.));

		for(std::size_t idx=0; idx<voronoi.size(); ++idx)
		{
			// voronoi vertex and its corresponding delaunay triangle
			const t_vec& voronoivert = voronoi[idx];
			const auto& thetriag = triags[idx];

			std::vector<const t_vec*> neighbourverts;
			for(std::size_t neighbourIdx : neighbours[idx])
			{
				const t_vec& neighbourvert = voronoi[neighbourIdx];
				neighbourverts.push_back(&neighbourvert);

				QLineF line{QPointF{voronoivert[0], voronoivert[1]}, QPointF{neighbourvert[0], neighbourvert[1]}};
				QGraphicsItem *item = m_scene->addLine(line, penVoronoi);
				m_voronoi.insert(item);
			}

			// not all triangle edges have neighbours -> there are unbound regions
			if(neighbourverts.size() < 3)
			{
				//const bool voronoivert_in_triag = pt_inside_hull<t_vec>(thetriag, voronoivert);

				// slopes of existing voronoi edges
				std::vector<t_real> slopes;
				for(const t_vec* vec : neighbourverts)
					slopes.push_back(line_angle(voronoivert, *vec));

				// iterate delaunay triangle vertices
				for(std::size_t idx1=0; idx1<thetriag.size(); ++idx1)
				{
					std::size_t idx2 = idx1+1;
					if(idx2 >= thetriag.size())
						idx2 = 0;

					t_vec vecMid = thetriag[idx1] + (thetriag[idx2] - thetriag[idx1]) * t_real{0.5};
					t_real angle = line_angle(voronoivert, vecMid);

					// if the slope angle doesn't exist yet, it leads to an unbound external region
					if(auto iterSlope = std::find_if(slopes.begin(), slopes.end(), [angle](t_real angle2) -> bool
						{ return m::angle_equals<t_real>(angle, angle2, g_eps, m::pi<t_real>); });
						iterSlope == slopes.end())
					{
						t_vec vecUnbound = (vecMid-voronoivert);
						t_real lengthUnbound = 1000. / m::norm(vecUnbound);
						t_vec vecOuter = voronoivert;

						// voronoi vertex on other side of edge?
						if(side_of_line<t_vec>(thetriag[idx1], thetriag[idx2], voronoivert) < 0.)
							vecOuter -= lengthUnbound*vecUnbound;
						else
							vecOuter += lengthUnbound*vecUnbound;

						QLineF line{QPointF{voronoivert[0], voronoivert[1]}, QPointF{vecOuter[0], vecOuter[1]}};
						QGraphicsItem *item = m_scene->addLine(line, penVoronoiUnbound);
						m_voronoi.insert(item);
					}
				}
			}
		}
	}


	if(m_calcdelaunay)
	{
		QPen penDelaunay;
		penDelaunay.setStyle(Qt::SolidLine);
		penDelaunay.setWidthF(1.);
		penDelaunay.setColor(QColor::fromRgbF(0.,0.,0.));

		// delaunay triangles
		for(const auto& thetriag : triags)
		{
			for(std::size_t idx1=0; idx1<thetriag.size(); ++idx1)
			{
				std::size_t idx2 = idx1+1;
				if(idx2 >= thetriag.size())
					idx2 = 0;

				QLineF line{QPointF{thetriag[idx1][0], thetriag[idx1][1]}, QPointF{thetriag[idx2][0], thetriag[idx2][1]}};
				QGraphicsItem *item = m_scene->addLine(line, penDelaunay);
				m_delaunay.insert(item);
			}
		}
	}
}

// ----------------------------------------------------------------------------




// ----------------------------------------------------------------------------

HullWnd::HullWnd(QWidget* pParent) : QMainWindow{pParent},
	m_scene{new QGraphicsScene{this}},
	m_view{new HullView{m_scene.get(), this}},
	m_statusLabel{std::make_shared<QLabel>(this)}
{
	// ------------------------------------------------------------------------
	// restore settings
	QSettings settings{this};

	if(settings.contains("wnd_geo"))
	{
		QByteArray arr{settings.value("wnd_geo").toByteArray()};
		this->restoreGeometry(arr);
	}
	if(settings.contains("wnd_state"))
	{
		QByteArray arr{settings.value("wnd_state").toByteArray()};
		this->restoreState(arr);
	}

	m_view->SetCalculateHull(
		settings.value("calc_hull", m_view->GetCalculateHull()).toBool());
	m_view->SetCalculateVoronoiVertices(
		settings.value("calc_voronoivertices", m_view->GetCalculateVoronoiVertices()).toBool());
	m_view->SetCalculateVoronoiRegions(
		settings.value("calc_voronoiregions", m_view->GetCalculateVoronoiRegions()).toBool());
	m_view->SetCalculateDelaunay(
		settings.value("calc_delaunay", m_view->GetCalculateDelaunay()).toBool());
	// ------------------------------------------------------------------------


	m_view->setRenderHints(QPainter::Antialiasing);

	setWindowTitle("Geo2D");
	setCentralWidget(m_view.get());

	QStatusBar *statusBar = new QStatusBar{this};
	statusBar->addPermanentWidget(m_statusLabel.get(), 1);
	setStatusBar(statusBar);


	// menu actions
	QAction *actionNew = new QAction{"New", this};
	connect(actionNew, &QAction::triggered, [this]()
		{ m_view->ClearVertices(); });

	QAction *actionLoad = new QAction{"Load...", this};
	connect(actionLoad, &QAction::triggered, [this]()
	{
		if(QString file = QFileDialog::getOpenFileName(this, "Load Data", "",
			"XML Files (*.xml);;All Files (* *.*)"); file!="")
		{
			std::ifstream ifstr(file.toStdString());
			if(!ifstr)
			{
				QMessageBox::critical(this, "Error", "File could not be opened for loading.");
				return;
			}

			ptree::ptree prop{};
			ptree::read_xml(ifstr, prop);

			std::size_t vertidx = 0;
			while(true)
			{
				std::ostringstream ostrVert;
				ostrVert << "geo2d.hull.vertex_" << vertidx;

				auto vertprop = prop.get_child_optional(ostrVert.str());
				if(!vertprop)
					break;

				auto vertx = vertprop->get_optional<t_real>("<xmlattr>.x");
				auto verty = vertprop->get_optional<t_real>("<xmlattr>.y");

				if(!vertx || !verty)
					break;

				m_view->AddVertex(QPointF{*vertx, *verty});

				++vertidx;
			}

			if(vertidx > 0)
				m_view->UpdateAll();
			else
				QMessageBox::warning(this, "Warning", "File contains no data.");
		}
	});

	QAction *actionSaveAs = new QAction{"Save as...", this};
	connect(actionSaveAs, &QAction::triggered, [this]()
	{
		if(QString file = QFileDialog::getSaveFileName(this, "Save Data", "",
			"XML Files (*.xml);;All Files (* *.*)"); file!="")
		{
			std::ofstream ofstr(file.toStdString());
			if(!ofstr)
			{
				QMessageBox::critical(this, "Error", "File could not be opened for saving.");
				return;
			}

			ptree::ptree prop{};

			std::size_t vertidx = 0;
			for(const Vertex* vertex : m_view->GetVertices())
			{
				QPointF vertexpos = vertex->scenePos();

				std::ostringstream ostrX, ostrY;
				ostrX << "geo2d.hull.vertex_" << vertidx << ".<xmlattr>.x";
				ostrY << "geo2d.hull.vertex_" << vertidx << ".<xmlattr>.y";

				prop.put<t_real>(ostrX.str(), vertexpos.x());
				prop.put<t_real>(ostrY.str(), vertexpos.y());

				++vertidx;
			}

			ptree::write_xml(ofstr, prop, ptree::xml_writer_make_settings('\t', 1, std::string{"utf-8"}));
		}
	});

	QAction *actionQuit = new QAction{"Exit", this};
	connect(actionQuit, &QAction::triggered, [this]()
		{ this->close(); });


	QAction *actionHull = new QAction{"Convex Hull", this};
	actionHull->setCheckable(true);
	actionHull->setChecked(m_view->GetCalculateHull());
	connect(actionHull, &QAction::toggled, [this](bool b)
		{ m_view->SetCalculateHull(b); });

	QAction *actionVoronoi = new QAction{"Voronoi Vertices", this};
	actionVoronoi->setCheckable(true);
	actionVoronoi->setChecked(m_view->GetCalculateVoronoiVertices());
	connect(actionVoronoi, &QAction::toggled, [this](bool b)
		{ m_view->SetCalculateVoronoiVertices(b); });

	QAction *actionVoronoiRegions = new QAction{"Voronoi Regions", this};
	actionVoronoiRegions->setCheckable(true);
	actionVoronoiRegions->setChecked(m_view->GetCalculateVoronoiRegions());
	connect(actionVoronoiRegions, &QAction::toggled, [this](bool b)
		{ m_view->SetCalculateVoronoiRegions(b); });

	QAction *actionDelaunay = new QAction{"Delaunay Triangulation", this};
	actionDelaunay->setCheckable(true);
	actionDelaunay->setChecked(m_view->GetCalculateDelaunay());
	connect(actionDelaunay, &QAction::toggled, [this](bool b)
		{ m_view->SetCalculateDelaunay(b); });


	QAction *actionHullQHull = new QAction{"QHull", this};
	actionHullQHull->setCheckable(true);
	actionHullQHull->setChecked(true);
	connect(actionHullQHull, &QAction::toggled, [this]()
		{ m_view->SetHullCalculationMethod(HullCalculationMethod::QHULL); });

	QAction *actionHullContour = new QAction{"Contour", this};
	actionHullContour->setCheckable(true);
	connect(actionHullContour, &QAction::toggled, [this]()
		{ m_view->SetHullCalculationMethod(HullCalculationMethod::CONTOUR); });

	QAction *actionHullInc = new QAction{"Incremental", this};
	actionHullInc->setCheckable(true);
	connect(actionHullInc, &QAction::toggled, [this]()
	{ m_view->SetHullCalculationMethod(HullCalculationMethod::ITERATIVE); });

	QAction *actionHullDivide = new QAction{"Divide && Conquer", this};
	actionHullDivide->setCheckable(true);
	connect(actionHullDivide, &QAction::toggled, [this]()
		{ m_view->SetHullCalculationMethod(HullCalculationMethod::RECURSIVE); });


	QAction *actionDelaunayQHull = new QAction{"QHull", this};
	actionDelaunayQHull->setCheckable(true);
	actionDelaunayQHull->setChecked(true);
	connect(actionDelaunayQHull, &QAction::toggled, [this]()
		{ m_view->SetDelaunayCalculationMethod(DelaunayCalculationMethod::QHULL); });

	QAction *actionDelaunayPara = new QAction{"Parabolic Trafo", this};
	actionDelaunayPara->setCheckable(true);
	connect(actionDelaunayPara, &QAction::toggled, [this]()
		{ m_view->SetDelaunayCalculationMethod(DelaunayCalculationMethod::PARABOLIC); });


	QActionGroup *groupHullBack = new QActionGroup{this};
	groupHullBack->addAction(actionHullQHull);
	groupHullBack->addAction(actionHullContour);
	groupHullBack->addAction(actionHullInc);
	groupHullBack->addAction(actionHullDivide);

	QActionGroup *groupDelaunayBack = new QActionGroup{this};
	groupDelaunayBack->addAction(actionDelaunayQHull);
	groupDelaunayBack->addAction(actionDelaunayPara);


	// menu
	QMenu *menuFile = new QMenu{"File", this};
	QMenu *menuCalc = new QMenu{"Calculate", this};
	QMenu *menuBack = new QMenu{"Backend", this};

	menuFile->addAction(actionNew);
	menuFile->addSeparator();
	menuFile->addAction(actionLoad);
	menuFile->addAction(actionSaveAs);
	menuFile->addSeparator();
	menuFile->addAction(actionQuit);

	menuCalc->addAction(actionHull);
	menuCalc->addSeparator();
	menuCalc->addAction(actionVoronoi);
	menuCalc->addAction(actionVoronoiRegions);
	menuCalc->addSeparator();
	menuCalc->addAction(actionDelaunay);

	menuBack->addSeparator()->setText("Convex Hull");
	menuBack->addAction(actionHullQHull);
	menuBack->addAction(actionHullContour);
	menuBack->addAction(actionHullInc);
	menuBack->addAction(actionHullDivide);
	menuBack->addSeparator()->setText("Delaunay");
	menuBack->addAction(actionDelaunayQHull);
	menuBack->addAction(actionDelaunayPara);


	// menu bar
	QMenuBar *menuBar = new QMenuBar{this};
	menuBar->addMenu(menuFile);
	menuBar->addMenu(menuCalc);
	menuBar->addMenu(menuBack);
	setMenuBar(menuBar);


	// connections
	connect(m_view.get(), &HullView::SignalMouseCoordinates, [this](double x, double y) -> void
	{
		SetStatusMessage(QString("x=%1, y=%2.").arg(x, 5).arg(y, 5));
	});


	SetStatusMessage("Ready.");
}


void HullWnd::SetStatusMessage(const QString& msg)
{
	m_statusLabel->setText(msg);
}


void HullWnd::closeEvent(QCloseEvent *e)
{
	// ------------------------------------------------------------------------
	// save settings
	QSettings settings{this};

	QByteArray geo{this->saveGeometry()}, state{this->saveState()};
	settings.setValue("wnd_geo", geo);
	settings.setValue("wnd_state", state);
	settings.setValue("calc_hull", m_view->GetCalculateHull());
	settings.setValue("calc_voronoivertices", m_view->GetCalculateVoronoiVertices());
	settings.setValue("calc_voronoiregions", m_view->GetCalculateVoronoiRegions());
	settings.setValue("calc_delaunay", m_view->GetCalculateDelaunay());
	// ------------------------------------------------------------------------

	QMainWindow::closeEvent(e);
}


HullWnd::~HullWnd()
{
}


// ----------------------------------------------------------------------------




// ----------------------------------------------------------------------------

static inline void set_locales()
{
	std::ios_base::sync_with_stdio(false);

	::setlocale(LC_ALL, "C");
	std::locale::global(std::locale("C"));
	QLocale::setDefault(QLocale::C);
}


int main(int argc, char** argv)
{
	try
	{
		auto app = std::make_unique<QApplication>(argc, argv);
		app->setOrganizationName("tw");
		app->setApplicationName("geo2d");
		set_locales();

		auto hullwnd = std::make_unique<HullWnd>();
		hullwnd->resize(1024, 768);
		hullwnd->show();

		return app->exec();
	}
	catch(const std::exception& ex)
	{
		std::cerr << ex.what() << std::endl;
	}

	return -1;
}
// ----------------------------------------------------------------------------

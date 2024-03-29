/**
 * convex hull test program
 * @author Tobias Weber (orcid: 0000-0002-7230-1932)
 * @date 15-Aug-2020
 * @license: see 'LICENSE' file
 */

#include "voro_gui.h"

#include <QApplication>
#include <QMenuBar>
#include <QLabel>
#include <QStatusBar>
#include <QMouseEvent>
#include <QFileDialog>
#include <QSvgGenerator>
#include <QMessageBox>
#include <QSettings>

#include <locale>
#include <memory>
#include <array>
#include <vector>
#include <iostream>

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>
namespace ptree = boost::property_tree;

#include "geo_algos.h"

using t_real = double;
using t_vec = m::vec<t_real, std::vector>;
using t_mat = m::mat<t_real, std::vector>;

const constexpr t_real g_eps = 1e-5;


//#define HULL_CHECK


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

HullScene::HullScene(QWidget* parent) : QGraphicsScene(parent), m_parent{parent}
{
}

HullScene::~HullScene()
{
}


void HullScene::AddVertex(const QPointF& pos)
{
	Vertex *vertex = new Vertex{pos};
	m_vertices.insert(vertex);
	addItem(vertex);
}



void HullScene::SetCalculateHull(bool b)
{
	m_calchull = b;
	UpdateHull();
}


void HullScene::SetCalculateVoronoiVertices(bool b)
{
	m_calcvoronoivertices = b;
	UpdateDelaunay();
}


void HullScene::SetCalculateVoronoiCircles(bool b)
{
	m_calcvoronoicircles = b;
	UpdateDelaunay();
}


void HullScene::SetCalculateVoronoiRegions(bool b)
{
	m_calcvoronoiregions = b;
	UpdateDelaunay();
}


void HullScene::SetCalculateDelaunay(bool b)
{
	m_calcdelaunay = b;
	UpdateDelaunay();
}


void HullScene::SetCalculateKruskal(bool b)
{
	m_calckruskal = b;
	UpdateDelaunay();
}


void HullScene::SetHullCalculationMethod(HullCalculationMethod m)
{
	m_hullcalculationmethod = m;
	UpdateHull();
}


void HullScene::SetDelaunayCalculationMethod(DelaunayCalculationMethod m)
{
	m_delaunaycalculationmethod = m;
	UpdateDelaunay();
}


void HullScene::SetSpanCalculationMethod(SpanCalculationMethod m)
{
	m_spancalculationmethod = m;
	UpdateDelaunay();
}


void HullScene::ClearVertices()
{
	for(Vertex* vertex : m_vertices)
	{
		removeItem(vertex);
		delete vertex;
	}
	m_vertices.clear();

	UpdateAll();
}


void HullScene::UpdateAll()
{
	UpdateDelaunay();
	UpdateHull();
}


void HullScene::UpdateHull()
{
	// remove previous hull
	for(QGraphicsItem* hullItem : m_hull)
	{
		removeItem(hullItem);
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
			std::tie(std::ignore, hull, std::ignore) = g::calc_delaunay<t_vec>(2, vertices, true);
			break;
		case HullCalculationMethod::CONTOUR:
			hull.emplace_back(g::calc_hull_contour<t_vec>(vertices, g_eps));
			break;
		case HullCalculationMethod::ITERATIVE:
			//hull.emplace_back(calc_hull_iterative<t_vec>(vertices));
			hull.emplace_back(g::calc_hull_iterative_bintree<t_vec>(vertices, g_eps));
			break;
		case HullCalculationMethod::RECURSIVE:
			hull.emplace_back(g::calc_hull_recursive<t_vec>(vertices, g_eps));
			break;
		default:
			QMessageBox::critical(m_parent, "Error", "Unknown hull calculation method.");
			break;
	}


#ifdef HULL_CHECK
	std::vector<t_vec> hullvertices;
	for(const auto& thetriag : hull)
		for(std::size_t idx1=0; idx1<thetriag.size(); ++idx1)
			hullvertices.emplace_back(thetriag[idx1]);
#endif

	// convex hull
	QPen penHull;
	penHull.setWidthF(3.);
	penHull.setStyle(Qt::SolidLine);
	penHull.setColor(QColor::fromRgbF(0.,0.,1.));

	for(const auto& thetriag : hull)
	{
		for(std::size_t idx1=0; idx1<thetriag.size(); ++idx1)
		{
			std::size_t idx2 = idx1+1;
			if(idx2 >= thetriag.size())
				idx2 = 0;
			if(idx1 == idx2)
				continue;

			#ifdef HULL_CHECK
			if(!all_points_on_same_side(thetriag[idx1], thetriag[idx2], hullvertices, g_eps))
				continue;
			#endif

			QLineF line{QPointF{thetriag[idx1][0], thetriag[idx1][1]}, QPointF{thetriag[idx2][0], thetriag[idx2][1]}};
			QGraphicsItem *item = addLine(line, penHull);
			m_hull.insert(item);
		}
	}
}


void HullScene::UpdateDelaunay()
{
	// remove previous triangulation
	for(QGraphicsItem* item : m_delaunay)
	{
		removeItem(item);
		delete item;
	}
	m_delaunay.clear();

	// remove previous voronoi vertices
	for(QGraphicsItem* item : m_voronoi)
	{
		removeItem(item);
		delete item;
	}
	m_voronoi.clear();


	if((!m_calcdelaunay && !m_calckruskal
		&& !m_calcvoronoivertices && !m_calcvoronoicircles && !m_calcvoronoiregions)
		|| m_vertices.size() < 4)
		return;


	// get vertices
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
			std::tie(voronoi, triags, neighbours) = g::calc_delaunay<t_vec>(2, vertices, false);
			break;
		case DelaunayCalculationMethod::ITERATIVE:
			std::tie(voronoi, triags, neighbours) = g::calc_delaunay_iterative<t_vec>(vertices, g_eps);
			break;
		case DelaunayCalculationMethod::PARABOLIC:
			std::tie(voronoi, triags, neighbours) = g::calc_delaunay_parabolic<t_vec>(vertices);
			break;
		default:
			QMessageBox::critical(m_parent, "Error", "Unknown Delaunay calculation method.");
			break;
	}


	const t_real itemRad = 7.;

	if(m_calcvoronoivertices || m_calcvoronoicircles)
	{
		QPen penVoronoi;
		penVoronoi.setStyle(Qt::SolidLine);
		penVoronoi.setWidthF(1.);
		penVoronoi.setColor(QColor::fromRgbF(1.,0.,0.));

		QPen penCircle;
		penCircle.setStyle(Qt::DotLine);
		penCircle.setWidthF(1.);
		penCircle.setColor(QColor::fromRgbF(1.,0.,0.));

		QBrush brushVoronoi;
		brushVoronoi.setStyle(Qt::SolidPattern);
		brushVoronoi.setColor(QColor::fromRgbF(1.,0.,0.));

		for(std::size_t idx=0; idx<voronoi.size(); ++idx)
		{
			const t_vec& voronoivert = voronoi[idx];
			QPointF voronoipt{voronoivert[0], voronoivert[1]};

			// voronoi vertices
			if(m_calcvoronoivertices)
			{
				QGraphicsItem *voronoiItem = addEllipse(
					voronoipt.x()-itemRad/2., voronoipt.y()-itemRad/2., itemRad, itemRad, penVoronoi, brushVoronoi);
				m_voronoi.insert(voronoiItem);
			}

			// circles
			if(m_calcvoronoicircles && idx < triags.size())
			{
				const auto& triag = triags[idx];
				if(triag.size() >= 3)
				{
					t_real rad = m::norm(voronoivert-triag[0]);

					QGraphicsItem *voronoiCircle = addEllipse(
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
		penVoronoi.setColor(QColor::fromRgbF(0.,0.,0.));

		QPen penVoronoiUnbound;
		penVoronoiUnbound.setStyle(Qt::DashLine);
		penVoronoiUnbound.setWidthF(1.);
		penVoronoiUnbound.setColor(QColor::fromRgbF(0.,0.,0.));

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
				QGraphicsItem *item = addLine(line, penVoronoi);
				m_voronoi.insert(item);
			}

			// not all triangle edges have neighbours -> there are unbound regions
			if(neighbourverts.size() < 3)
			{
				//const bool voronoivert_in_triag = pt_inside_hull<t_vec>(thetriag, voronoivert);

				// slopes of existing voronoi edges
				std::vector<t_real> slopes;
				for(const t_vec* vec : neighbourverts)
					slopes.push_back(g::line_angle(voronoivert, *vec));

				// iterate delaunay triangle vertices
				for(std::size_t idx1=0; idx1<thetriag.size(); ++idx1)
				{
					std::size_t idx2 = idx1+1;
					if(idx2 >= thetriag.size())
						idx2 = 0;

					t_vec vecMid = thetriag[idx1] + (thetriag[idx2] - thetriag[idx1]) * t_real{0.5};
					t_real angle = g::line_angle(voronoivert, vecMid);

					// if the slope angle doesn't exist yet, it leads to an unbound external region
					if(auto iterSlope = std::find_if(slopes.begin(), slopes.end(), [angle](t_real angle2) -> bool
					{ return m::angle_equals<t_real>(angle, angle2, g_eps, m::pi<t_real>); });
					iterSlope == slopes.end())
					{
						t_vec vecUnbound = (vecMid-voronoivert);
						t_real lengthUnbound = 1000. / m::norm(vecUnbound);
						t_vec vecOuter = voronoivert;

						// voronoi vertex on other side of edge?
						if(g::side_of_line<t_vec>(thetriag[idx1], thetriag[idx2], voronoivert) < 0.)
							vecOuter -= lengthUnbound*vecUnbound;
						else
							vecOuter += lengthUnbound*vecUnbound;

						QLineF line{QPointF{voronoivert[0], voronoivert[1]}, QPointF{vecOuter[0], vecOuter[1]}};
						QGraphicsItem *item = addLine(line, penVoronoiUnbound);
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
		penDelaunay.setColor(QColor::fromRgbF(0.,0.,1.));

		// delaunay triangles
		for(const auto& thetriag : triags)
		{
			for(std::size_t idx1=0; idx1<thetriag.size(); ++idx1)
			{
				std::size_t idx2 = idx1+1;
				if(idx2 >= thetriag.size())
					idx2 = 0;

				QLineF line{QPointF{thetriag[idx1][0], thetriag[idx1][1]}, QPointF{thetriag[idx2][0], thetriag[idx2][1]}};
				QGraphicsItem *item = addLine(line, penDelaunay);
				m_delaunay.insert(item);
			}
		}
	}


	if(m_calckruskal)
	{
		QPen penKruskal;
		penKruskal.setStyle(Qt::SolidLine);
		penKruskal.setWidthF(2.);
		penKruskal.setColor(QColor::fromRgbF(0.,0.7,0.));

		auto edges = g::get_edges(vertices, triags, g_eps);
		std::vector<std::pair<std::size_t, std::size_t>> span;

		switch(m_spancalculationmethod)
		{
			case SpanCalculationMethod::KRUSKAL:
				span = g::calc_min_spantree<t_vec>(vertices, edges);
				break;
			case SpanCalculationMethod::BOOST:
				span = g::calc_min_spantree_boost<t_vec>(vertices);
				break;
			default:
				QMessageBox::critical(m_parent, "Error", "Unknown span tree calculation method.");
				break;
		}

		for(const auto& spanedge : span)
		{
			const t_vec& vert1 = vertices[spanedge.first];
			const t_vec& vert2 = vertices[spanedge.second];

			QLineF line{QPointF{vert1[0], vert1[1]}, QPointF{vert2[0], vert2[1]}};
			QGraphicsItem *item = addLine(line, penKruskal);
			m_delaunay.insert(item);
		}
	}
}


// ----------------------------------------------------------------------------




// ----------------------------------------------------------------------------

HullView::HullView(HullScene *scene, QWidget *parent) : QGraphicsView(scene, parent),
	m_scene{scene}
{
	setHorizontalScrollBarPolicy(Qt::ScrollBarAlwaysOn);
	setVerticalScrollBarPolicy(Qt::ScrollBarAlwaysOn);

	setInteractive(true);
	setMouseTracking(true);

	//scale(1., -1.);

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
	for(const Vertex* vertex : m_scene->GetVertices())
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



void HullView::mousePressEvent(QMouseEvent *evt)
{
	QPoint posVP = evt->pos();
	QPointF posScene = mapToScene(posVP);

	QList<QGraphicsItem*> items = this->items(posVP);
	QGraphicsItem* item = nullptr;
	bool item_is_vertex = false;

	auto& verts = m_scene->GetVertices();

	for(int itemidx=0; itemidx<items.size(); ++itemidx)
	{
		item = items[itemidx];
		item_is_vertex = verts.find(static_cast<Vertex*>(item)) != verts.end();
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
			m_scene->AddVertex(posScene);
			m_dragging = true;
			m_scene->UpdateAll();
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
			verts.erase(static_cast<Vertex*>(item));
			delete item;
			m_scene->UpdateAll();
		}
	}

	QGraphicsView::mousePressEvent(evt);
}


void HullView::mouseReleaseEvent(QMouseEvent *evt)
{
	if(evt->button() == Qt::LeftButton)
		m_dragging = false;

	m_scene->UpdateAll();
	QGraphicsView::mouseReleaseEvent(evt);
}


void HullView::mouseMoveEvent(QMouseEvent *evt)
{
	QGraphicsView::mouseMoveEvent(evt);

	if(m_dragging)
	{
		QResizeEvent evt{size(), size()};
		resizeEvent(&evt);
		m_scene->UpdateAll();
	}

	QPoint posVP = evt->pos();
	QPointF posScene = mapToScene(posVP);
	emit SignalMouseCoordinates(posScene.x(), posScene.y());
}

// ----------------------------------------------------------------------------




// ----------------------------------------------------------------------------

HullWnd::HullWnd(QWidget* pParent) : QMainWindow{pParent},
	m_scene{new HullScene{this}},
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
	else
	{
		resize(1024, 768);
	}
	if(settings.contains("wnd_state"))
	{
		QByteArray arr{settings.value("wnd_state").toByteArray()};
		this->restoreState(arr);
	}

	m_scene->SetCalculateHull(
		settings.value("calc_hull", m_scene->GetCalculateHull()).toBool());
	m_scene->SetCalculateVoronoiVertices(
		settings.value("calc_voronoivertices", m_scene->GetCalculateVoronoiVertices()).toBool());
	m_scene->SetCalculateVoronoiCircles(
		settings.value("calc_voronoicircles", m_scene->GetCalculateVoronoiCircles()).toBool());
	m_scene->SetCalculateVoronoiRegions(
		settings.value("calc_voronoiregions", m_scene->GetCalculateVoronoiRegions()).toBool());
	m_scene->SetCalculateDelaunay(
		settings.value("calc_delaunay", m_scene->GetCalculateDelaunay()).toBool());
	m_scene->SetCalculateKruskal(
		settings.value("calc_kruskal", m_scene->GetCalculateKruskal()).toBool());
	// ------------------------------------------------------------------------


	m_view->setRenderHints(QPainter::Antialiasing);

	setWindowTitle("Voronoi Diagrams");
	setCentralWidget(m_view.get());

	QStatusBar *statusBar = new QStatusBar{this};
	statusBar->addPermanentWidget(m_statusLabel.get(), 1);
	setStatusBar(statusBar);


	// menu actions
	QAction *actionNew = new QAction{"New", this};
	connect(actionNew, &QAction::triggered, [this]()
		{ m_scene->ClearVertices(); });

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

			m_scene->ClearVertices();

			ptree::ptree prop{};
			ptree::read_xml(ifstr, prop);

			std::size_t vertidx = 0;
			while(true)
			{
				std::ostringstream ostrVert;
				ostrVert << "voro2d.vertices." << vertidx;

				auto vertprop = prop.get_child_optional(ostrVert.str());
				if(!vertprop)
					break;

				auto vertx = vertprop->get_optional<t_real>("<xmlattr>.x");
				auto verty = vertprop->get_optional<t_real>("<xmlattr>.y");

				if(!vertx || !verty)
					break;

				m_scene->AddVertex(QPointF{*vertx, *verty});

				++vertidx;
			}

			if(vertidx > 0)
				m_scene->UpdateAll();
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
			for(const Vertex* vertex : m_scene->GetVertices())
			{
				QPointF vertexpos = vertex->scenePos();

				std::ostringstream ostrX, ostrY;
				ostrX << "voro2d.vertices." << vertidx << ".<xmlattr>.x";
				ostrY << "voro2d.vertices." << vertidx << ".<xmlattr>.y";

				prop.put<t_real>(ostrX.str(), vertexpos.x());
				prop.put<t_real>(ostrY.str(), vertexpos.y());

				++vertidx;
			}

			ptree::write_xml(ofstr, prop, ptree::xml_writer_make_settings('\t', 1, std::string{"utf-8"}));
		}
	});

	QAction *actionExportSvg = new QAction{"Export SVG...", this};
	connect(actionExportSvg, &QAction::triggered, [this]()
	{
		if(QString file = QFileDialog::getSaveFileName(this, "Export SVG", "",
			"SVG Files (*.svg);;All Files (* *.*)"); file!="")
		{
			QSvgGenerator svggen;
			svggen.setSize(QSize{width(), height()});
			svggen.setFileName(file);

			QPainter paint(&svggen);
			m_scene->render(&paint);
		}
	});

	QAction *actionQuit = new QAction{"Exit", this};
	connect(actionQuit, &QAction::triggered, [this]()
		{ this->close(); });


	QAction *actionHull = new QAction{"Convex Hull", this};
	actionHull->setCheckable(true);
	actionHull->setChecked(m_scene->GetCalculateHull());
	connect(actionHull, &QAction::toggled, [this](bool b)
		{ m_scene->SetCalculateHull(b); });

	QAction *actionVoronoi = new QAction{"Voronoi Vertices", this};
	actionVoronoi->setCheckable(true);
	actionVoronoi->setChecked(m_scene->GetCalculateVoronoiVertices());
	connect(actionVoronoi, &QAction::toggled, [this](bool b)
		{ m_scene->SetCalculateVoronoiVertices(b); });

	QAction *actionVoronoiCircles = new QAction{"Voronoi Circles", this};
	actionVoronoiCircles->setCheckable(true);
	actionVoronoiCircles->setChecked(m_scene->GetCalculateVoronoiCircles());
	connect(actionVoronoiCircles, &QAction::toggled, [this](bool b)
		{ m_scene->SetCalculateVoronoiCircles(b); });

	QAction *actionVoronoiRegions = new QAction{"Voronoi Regions", this};
	actionVoronoiRegions->setCheckable(true);
	actionVoronoiRegions->setChecked(m_scene->GetCalculateVoronoiRegions());
	connect(actionVoronoiRegions, &QAction::toggled, [this](bool b)
		{ m_scene->SetCalculateVoronoiRegions(b); });

	QAction *actionDelaunay = new QAction{"Delaunay Triangulation", this};
	actionDelaunay->setCheckable(true);
	actionDelaunay->setChecked(m_scene->GetCalculateDelaunay());
	connect(actionDelaunay, &QAction::toggled, [this](bool b)
		{ m_scene->SetCalculateDelaunay(b); });

	QAction *actionSpanTree = new QAction{"Minimum Spanning Tree", this};
	actionSpanTree->setCheckable(true);
	actionSpanTree->setChecked(m_scene->GetCalculateKruskal());
	connect(actionSpanTree, &QAction::toggled, [this](bool b)
		{ m_scene->SetCalculateKruskal(b); });


	QAction *actionHullQHull = new QAction{"QHull", this};
	actionHullQHull->setCheckable(true);
	actionHullQHull->setChecked(true);
	connect(actionHullQHull, &QAction::toggled, [this]()
		{ m_scene->SetHullCalculationMethod(HullCalculationMethod::QHULL); });

	QAction *actionHullContour = new QAction{"Contour", this};
	actionHullContour->setCheckable(true);
	connect(actionHullContour, &QAction::toggled, [this]()
		{ m_scene->SetHullCalculationMethod(HullCalculationMethod::CONTOUR); });

	QAction *actionHullInc = new QAction{"Incremental", this};
	actionHullInc->setCheckable(true);
	connect(actionHullInc, &QAction::toggled, [this]()
		{ m_scene->SetHullCalculationMethod(HullCalculationMethod::ITERATIVE); });

	QAction *actionHullDivide = new QAction{"Divide && Conquer", this};
	actionHullDivide->setCheckable(true);
	connect(actionHullDivide, &QAction::toggled, [this]()
		{ m_scene->SetHullCalculationMethod(HullCalculationMethod::RECURSIVE); });


	QAction *actionDelaunayQHull = new QAction{"QHull", this};
	actionDelaunayQHull->setCheckable(true);
	actionDelaunayQHull->setChecked(true);
	connect(actionDelaunayQHull, &QAction::toggled, [this]()
		{ m_scene->SetDelaunayCalculationMethod(DelaunayCalculationMethod::QHULL); });

	QAction *actionDelaunayInc = new QAction{"Incremental", this};
	actionDelaunayInc->setCheckable(true);
	connect(actionDelaunayInc, &QAction::toggled, [this]()
		{ m_scene->SetDelaunayCalculationMethod(DelaunayCalculationMethod::ITERATIVE); });

	QAction *actionDelaunayPara = new QAction{"Parabolic Trafo", this};
	actionDelaunayPara->setCheckable(true);
	connect(actionDelaunayPara, &QAction::toggled, [this]()
		{ m_scene->SetDelaunayCalculationMethod(DelaunayCalculationMethod::PARABOLIC); });


	QAction *actionSpanKruskal = new QAction{"Kruskal", this};
	actionSpanKruskal->setCheckable(true);
	actionSpanKruskal->setChecked(true);
	connect(actionSpanKruskal, &QAction::toggled, [this]()
		{ m_scene->SetSpanCalculationMethod(SpanCalculationMethod::KRUSKAL); });

	QAction *actionSpanBoost = new QAction{"Kruskal via Boost.Graph", this};
	actionSpanBoost->setCheckable(true);
	connect(actionSpanBoost, &QAction::toggled, [this]()
		{ m_scene->SetSpanCalculationMethod(SpanCalculationMethod::BOOST); });


	QActionGroup *groupHullBack = new QActionGroup{this};
	groupHullBack->addAction(actionHullQHull);
	groupHullBack->addAction(actionHullContour);
	groupHullBack->addAction(actionHullInc);
	groupHullBack->addAction(actionHullDivide);

	QActionGroup *groupDelaunayBack = new QActionGroup{this};
	groupDelaunayBack->addAction(actionDelaunayQHull);
	groupDelaunayBack->addAction(actionDelaunayInc);
	groupDelaunayBack->addAction(actionDelaunayPara);

	QActionGroup *groupSpanBack = new QActionGroup{this};
	groupSpanBack->addAction(actionSpanKruskal);
	groupSpanBack->addAction(actionSpanBoost);


	// menu
	QMenu *menuFile = new QMenu{"File", this};
	QMenu *menuCalc = new QMenu{"Calculate", this};
	QMenu *menuBack = new QMenu{"Backends", this};

	menuFile->addAction(actionNew);
	menuFile->addSeparator();
	menuFile->addAction(actionLoad);
	menuFile->addAction(actionSaveAs);
	menuFile->addSeparator();
	menuFile->addAction(actionExportSvg);
	menuFile->addSeparator();
	menuFile->addAction(actionQuit);

	menuCalc->addAction(actionHull);
	menuCalc->addSeparator();
	menuCalc->addAction(actionVoronoi);
	menuCalc->addAction(actionVoronoiCircles);
	menuCalc->addAction(actionVoronoiRegions);
	menuCalc->addSeparator();
	menuCalc->addAction(actionDelaunay);
	menuCalc->addAction(actionSpanTree);


	QMenu *menuBackHull = new QMenu{"Convex Hull", this};
	menuBackHull->addAction(actionHullQHull);
	menuBackHull->addAction(actionHullContour);
	menuBackHull->addAction(actionHullInc);
	menuBackHull->addAction(actionHullDivide);

	QMenu *menuDelaunay = new QMenu{"Delaunay Triangulation", this};
	menuDelaunay->addAction(actionDelaunayQHull);
	menuDelaunay->addAction(actionDelaunayInc);
	menuDelaunay->addAction(actionDelaunayPara);

	QMenu *menuSpan = new QMenu{"Minimum Spanning Tree", this};
	menuSpan->addAction(actionSpanKruskal);
	menuSpan->addAction(actionSpanBoost);

	menuBack->addMenu(menuBackHull);
	menuBack->addMenu(menuDelaunay);
	menuBack->addMenu(menuSpan);


	// menu bar
	QMenuBar *menuBar = new QMenuBar{this};
	menuBar->setNativeMenuBar(false);
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
	settings.setValue("calc_hull", m_scene->GetCalculateHull());
	settings.setValue("calc_voronoivertices", m_scene->GetCalculateVoronoiVertices());
	settings.setValue("calc_voronoicircles", m_scene->GetCalculateVoronoiCircles());
	settings.setValue("calc_voronoiregions", m_scene->GetCalculateVoronoiRegions());
	settings.setValue("calc_delaunay", m_scene->GetCalculateDelaunay());
	settings.setValue("calc_kruskal", m_scene->GetCalculateKruskal());
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
		app->setApplicationName("voro2d");
		set_locales();

		auto hullwnd = std::make_unique<HullWnd>();
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

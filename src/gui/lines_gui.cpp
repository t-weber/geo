/**
 * line intersection test program
 * @author Tobias Weber (orcid: 0000-0002-7230-1932)
 * @date 11-Nov-2020
 * @license: see 'LICENSE' file
 */

#include "lines_gui.h"

#include <QApplication>
#include <QMenuBar>
#include <QLabel>
#include <QStatusBar>
#include <QMouseEvent>
#include <QFileDialog>
#include <QSvgGenerator>
#include <QMessageBox>
#include <QSettings>
#include <QProgressDialog>

#include <locale>
#include <memory>
#include <array>
#include <vector>
#include <unordered_map>
#include <thread>
#include <mutex>
#include <future>
#include <iostream>

#include <boost/asio.hpp>
namespace asio = boost::asio;

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>
namespace ptree = boost::property_tree;


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

LinesScene::LinesScene(QWidget *parent) : QGraphicsScene(parent), m_parent{parent}
{
	ClearVertices();
}


LinesScene::~LinesScene()
{
	if(m_elem_voro)
		delete m_elem_voro;
}


void LinesScene::CreateVoroImage(int width, int height)
{
	// only delete old image and create a new one if the sizes have changed
	if(m_elem_voro && m_elem_voro->width()!=width && m_elem_voro->height()!=height)
	{
		delete m_elem_voro;
		m_elem_voro = nullptr;
	}

	if(!m_elem_voro)
	{
		m_elem_voro = new QImage(width, height, QImage::Format_RGB32);
		m_elem_voro->fill(/*backgroundBrush().color()*/ QColor::fromRgbF(0.95, 0.95, 0.95, 1.));
	}
}


void LinesScene::AddVertex(const QPointF& pos)
{
	Vertex *vertex = new Vertex{pos};
	m_elems_vertices.push_back(vertex);
	addItem(vertex);
}


void LinesScene::ClearVertices()
{
	for(Vertex* vertex : m_elems_vertices)
	{
		removeItem(vertex);
		delete vertex;
	}
	m_elems_vertices.clear();

	setBackgroundBrush(QBrush{QColor::fromRgbF(0.95, 0.95, 0.95, 1.)});
	if(m_elem_voro)
		m_elem_voro->fill(backgroundBrush().color());
	UpdateAll();
}


void LinesScene::SetIntersectionCalculationMethod(IntersectionCalculationMethod m)
{
	m_intersectioncalculationmethod = m;
	UpdateIntersections();
}


void LinesScene::SetVoronoiCalculationMethod(VoronoiCalculationMethod m)
{
	m_voronoicalculationmethod = m;
	UpdateVoro();
}


void LinesScene::UpdateAll()
{
	UpdateLines();
	UpdateIntersections();
	UpdateVoro();
}


void LinesScene::UpdateLines()
{
	// remove previous lines
	for(QGraphicsItem* item : m_elems_lines)
	{
		removeItem(item);
		delete item;
	}
	m_elems_lines.clear();


	// get new lines
	m_lines.clear();
	if(m_elems_vertices.size() < 2)
		return;
	m_lines.reserve(m_elems_vertices.size()/2);

	for(std::size_t i=0; i<m_elems_vertices.size()-1; i+=2)
	{
		const Vertex* _vert1 = m_elems_vertices[i];
		const Vertex* _vert2 = m_elems_vertices[i+1];

		t_vec vert1 = m::create<t_vec>({_vert1->x(), _vert1->y()});
		t_vec vert2 = m::create<t_vec>({_vert2->x(), _vert2->y()});

		// ensure that first vertex is on the left-hand side
		//if(vert1[0] > vert2[0])
		//	std::swap(vert1, vert2);

		m_lines.emplace_back(std::make_pair(vert1, vert2));
	}


	QPen penEdge;
	penEdge.setStyle(Qt::SolidLine);
	penEdge.setWidthF(2.);
	penEdge.setColor(QColor::fromRgbF(0., 0., 1.));


	for(const auto& line : m_lines)
	{
		const t_vec& vertex1 = line.first;
		const t_vec& vertex2 = line.second;

		QLineF qline{QPointF{vertex1[0], vertex1[1]}, QPointF{vertex2[0], vertex2[1]}};
		QGraphicsItem *item = addLine(qline, penEdge);
		m_elems_lines.push_back(item);
	}
}


void LinesScene::UpdateIntersections()
{
	// remove previous intersection points
	for(QGraphicsItem* item : m_elems_inters)
	{
		removeItem(item);
		delete item;
	}
	m_elems_inters.clear();


	std::vector<std::tuple<std::size_t, std::size_t, t_vec>> intersections;

	switch(m_intersectioncalculationmethod)
	{
		case IntersectionCalculationMethod::DIRECT:
			intersections = g::intersect_ineff<t_vec, std::pair<t_vec, t_vec>>(m_lines, g_eps);
			break;
		case IntersectionCalculationMethod::SWEEP:
			intersections = g::intersect_sweep<t_vec, std::pair<t_vec, t_vec>>(m_lines, g_eps);
			break;
		default:
			QMessageBox::critical(m_parent, "Error", "Unknown intersection calculation method.");
			break;
	};


	QPen pen;
	pen.setStyle(Qt::SolidLine);
	pen.setWidthF(1.);
	pen.setColor(QColor::fromRgbF(0., 0.25, 0.));

	QBrush brush;
	brush.setStyle(Qt::SolidPattern);
	brush.setColor(QColor::fromRgbF(0., 0.75, 0.));

	for(const auto& intersection : intersections)
	{
		const t_vec& inters = std::get<2>(intersection);

		const t_real width = 14.;
		QRectF rect{inters[0]-width/2, inters[1]-width/2, width, width};
		QGraphicsItem *item = addEllipse(rect, pen, brush);
		m_elems_inters.push_back(item);
	}
}


void LinesScene::UpdateVoroImage(const QTransform& trafoSceneToVP)
{
	QTransform trafoVPToScene = trafoSceneToVP.inverted();

	if(!m_elem_voro)
		return;

	unsigned int num_threads = std::thread::hardware_concurrency();
	if(num_threads > 8)
		num_threads = 8;
	asio::thread_pool tp{num_threads};

	std::vector<std::shared_ptr<std::packaged_task<void()>>> packages;
	std::mutex mtx;

	const int width = m_elem_voro->width();
	const int height = m_elem_voro->height();
	std::unordered_map<std::size_t, QColor> linecolours;

	QProgressDialog progdlg(m_parent);
	progdlg.setWindowModality(Qt::WindowModal);
	progdlg.setMinimum(0);
	progdlg.setMaximum(height);
	QString msg = QString("Calculating Voronoi regions in %1 threads...").arg(num_threads);
	progdlg.setLabel(new QLabel(msg));

	for(int y=0; y<height; ++y)
	{
		auto package = std::make_shared<std::packaged_task<void()>>(
			[this, y, width, &linecolours, &mtx, &trafoVPToScene]() -> void
			{
				for(int x=0; x<width; ++x)
				{
					t_real scenex, sceney;
					trafoVPToScene.map(x, y, &scenex, &sceney);

					t_vec pt = m::create<t_vec>({scenex, sceney});
					std::size_t lineidx = GetClosestLineIdx(pt);

					// get colour for voronoi region
					QColor col{0xff, 0xff, 0xff, 0xff};

					std::lock_guard<std::mutex> _lck(mtx);
					auto iter = linecolours.find(lineidx);
					if(iter != linecolours.end())
					{
						col = iter->second;
					}
					else
					{
						col.setRgb(g::get_rand<int>(0,0xff), g::get_rand<int>(0,0xff), g::get_rand<int>(0,0xff));
						linecolours.insert(std::make_pair(lineidx, col));
					}

					m_elem_voro->setPixelColor(x, y, col);
				}
			});

		packages.push_back(package);
		asio::post(tp, [package]() -> void { if(package) (*package)(); });
	}

	for(int y=0; y<height; ++y)
	{
		if(progdlg.wasCanceled())
			break;
		progdlg.setValue(y);

		if(packages[y])
			packages[y]->get_future().get();
	}

	tp.join();
	progdlg.setValue(m_elem_voro->height());

	setBackgroundBrush(*m_elem_voro);
}


std::size_t LinesScene::GetClosestLineIdx(const t_vec& pt) const
{
	t_real mindist = std::numeric_limits<t_real>::max();
	std::size_t minidx = 0;

	for(std::size_t idx=0; idx<m_lines.size(); ++idx)
	{
		const auto& line = m_lines[idx];

		t_real dist = dist_pt_line(pt, line.first, line.second, false);
		if(dist < mindist)
		{
			mindist = dist;
			minidx = idx;
		}
	}

	return minidx;
}


void LinesScene::UpdateVoro()
{
	// remove previous voronoi diagram
	for(QGraphicsItem* item : m_elems_voro)
	{
		removeItem(item);
		delete item;
	}
	m_elems_voro.clear();

	if(!m_calcvoro)
		return;


	using t_line = std::pair<t_vec, t_vec>;

	// get vertices and bisectors
	std::vector<t_vec> vertices{};
	std::vector<t_line> linear_edges{}, linear_inf_edges{};
	std::vector<std::vector<t_vec>> all_parabolic_edges{};
	std::vector<t_line> linear_helper_edges{};

	switch(m_voronoicalculationmethod)
	{
		case VoronoiCalculationMethod::BOOST:
			std::tie(vertices, linear_edges, linear_inf_edges, all_parabolic_edges, linear_helper_edges, m_vorograph)
				= g::calc_voro<t_vec, t_line, decltype(m_vorograph)>(m_lines);
			break;
		case VoronoiCalculationMethod::CGAL:
			std::tie(vertices, linear_edges, linear_inf_edges, all_parabolic_edges, linear_helper_edges, m_vorograph)
				= g::calc_voro_cgal<t_vec, t_line, decltype(m_vorograph)>(m_lines);
			break;
	}

	// linear finite voronoi edges
	QPen penLinEdge;
	penLinEdge.setStyle(Qt::SolidLine);
	penLinEdge.setWidthF(1.);
	penLinEdge.setColor(QColor::fromRgbF(0.,0.,0.));

	for(const auto& linear_edge : linear_edges)
	{
		QLineF line{
			QPointF{std::get<0>(linear_edge)[0], std::get<0>(linear_edge)[1]},
			QPointF{std::get<1>(linear_edge)[0], std::get<1>(linear_edge)[1]} };
		QGraphicsItem *item = addLine(line, penLinEdge);
		m_elems_voro.push_back(item);
	}

	// linear infinite voronoi edges
	QPen penLinEdgeInf;
	penLinEdgeInf.setStyle(Qt::DashLine);
	penLinEdgeInf.setWidthF(1.);
	penLinEdgeInf.setColor(QColor::fromRgbF(0.,0.,0.));

	for(const auto& linear_edge : linear_inf_edges)
	{
		QLineF line{
			QPointF{std::get<0>(linear_edge)[0], std::get<0>(linear_edge)[1]},
			QPointF{std::get<1>(linear_edge)[0], std::get<1>(linear_edge)[1]} };
			QGraphicsItem *item = addLine(line, penLinEdgeInf);
			m_elems_voro.push_back(item);
	}

	// linear helper edges
	if(m_calcvorohelpers)
	{
		QPen penLinHelperEdge;
		penLinHelperEdge.setStyle(Qt::DotLine);
		penLinHelperEdge.setWidthF(1.);
		penLinHelperEdge.setColor(QColor::fromRgbF(1., 0., 0.));

		for(const auto& linear_edge : linear_helper_edges)
		{
			QLineF line{
				QPointF{std::get<0>(linear_edge)[0], std::get<0>(linear_edge)[1]},
				QPointF{std::get<1>(linear_edge)[0], std::get<1>(linear_edge)[1]} };
				QGraphicsItem *item = addLine(line, penLinHelperEdge);
				m_elems_voro.push_back(item);
		}
	}

	// parabolic voronoi edges
	QPen penParaEdge = penLinEdge;

	for(const auto& parabolic_edges : all_parabolic_edges)
	{
		QPolygonF poly;
		poly.reserve(parabolic_edges.size());
		for(const auto& parabolic_edge : parabolic_edges)
			poly << QPointF{parabolic_edge[0], parabolic_edge[1]};

		QPainterPath path;
		path.addPolygon(poly);

		QGraphicsItem *item = addPath(path, penParaEdge);
		m_elems_voro.push_back(item);
	}

	// voronoi vertices
	QPen penVertex;
	penVertex.setStyle(Qt::SolidLine);
	penVertex.setWidthF(1.);
	penVertex.setColor(QColor::fromRgbF(0.25, 0., 0.));

	QBrush brushVertex;
	brushVertex.setStyle(Qt::SolidPattern);
	brushVertex.setColor(QColor::fromRgbF(0.75, 0., 0.));

	for(const auto& vertex : vertices)
	{
		const t_real width = 8.;
		QRectF rect{vertex[0]-width/2, vertex[1]-width/2, width, width};
		QGraphicsItem *item = addEllipse(rect, penVertex, brushVertex);
		m_elems_voro.push_back(item);
	}
}

// ----------------------------------------------------------------------------



// ----------------------------------------------------------------------------

LinesView::LinesView(LinesScene *scene, QWidget *parent) : QGraphicsView(scene, parent),
	m_scene{scene}
{
	setHorizontalScrollBarPolicy(Qt::ScrollBarAlwaysOn);
	setVerticalScrollBarPolicy(Qt::ScrollBarAlwaysOn);

	setInteractive(true);
	setMouseTracking(true);

	//scale(1., -1.);
}


LinesView::~LinesView()
{
}


void LinesView::resizeEvent(QResizeEvent *evt)
{
	int widthView = evt->size().width();
	int heightView = evt->size().height();

	QPointF pt1{mapToScene(QPoint{0,0})};
	QPointF pt2{mapToScene(QPoint{widthView, heightView})};

	// include bounds given by vertices
	const double padding = 16;
	for(const Vertex* vertex : m_scene->GetVertexElems())
	{
		QPointF vertexpos = vertex->scenePos();

		if(vertexpos.x() < pt1.x())
			pt1.setX(vertexpos.x() - padding);
		if(vertexpos.x() > pt2.x())
			pt2.setX(vertexpos.x() + padding);
		if(vertexpos.y() < pt1.y())
			pt1.setY(vertexpos.y() - padding);
		if(vertexpos.y() > pt2.y())
			pt2.setY(vertexpos.y() + padding);
	}
	setSceneRect(QRectF{pt1, pt2});

	m_scene->CreateVoroImage(widthView, heightView);
}



void LinesView::mousePressEvent(QMouseEvent *evt)
{
	QPoint posVP = evt->pos();
	QPointF posScene = mapToScene(posVP);

	QList<QGraphicsItem*> items = this->items(posVP);
	QGraphicsItem* item = nullptr;
	bool item_is_vertex = false;
	auto &verts = m_scene->GetVertexElems();

	for(int itemidx=0; itemidx<items.size(); ++itemidx)
	{
		item = items[itemidx];
		auto iter = std::find(verts.begin(), verts.end(), static_cast<Vertex*>(item));
		item_is_vertex = (iter != verts.end());
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
			auto iter = std::find(verts.begin(), verts.end(), static_cast<Vertex*>(item));

			std::size_t idx = iter - verts.begin();
			if(iter != verts.end())
				iter = verts.erase(iter);
			delete item;

			// move remaining vertex of line to the end
			std::size_t otheridx = (idx % 2 == 0 ? idx : idx-1);
			if(otheridx < verts.size())
			{
				Vertex* vert = verts[otheridx];
				verts.erase(verts.begin()+otheridx);
				verts.push_back(vert);
			}

			m_scene->UpdateAll();
		}
	}

	QGraphicsView::mousePressEvent(evt);
}


void LinesView::mouseReleaseEvent(QMouseEvent *evt)
{
	if(evt->button() == Qt::LeftButton)
		m_dragging = false;

	m_scene->UpdateAll();
	QGraphicsView::mouseReleaseEvent(evt);
}


void LinesView::mouseMoveEvent(QMouseEvent *evt)
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
	emit SignalMouseCoordinates(posScene.x(), posScene.y(), posVP.x(), posVP.y());
}


void LinesView::drawBackground(QPainter* painter, const QRectF& rect)
{
	//QGraphicsView::drawBackground(painter, rect.intersected(sceneRect()));
	//QGraphicsView::drawBackground(painter, rect);

	// hack, because the background brush is drawn with respect to scene (0,0), not vp (0,0)
	// TODO: handle scene-viewport trafos other than translations
	if(m_scene->GetVoroImage())
		painter->drawImage(mapToScene(QPoint(0,0)), *m_scene->GetVoroImage()/*, mapFromScene(rect)*/);
	else
		QGraphicsView::drawBackground(painter, rect);
}
// ----------------------------------------------------------------------------



// ----------------------------------------------------------------------------

LinesWnd::LinesWnd(QWidget* pParent) : QMainWindow{pParent},
	m_scene{new LinesScene{this}},
	m_view{new LinesView{m_scene.get(), this}},
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
	// ------------------------------------------------------------------------


	m_view->setRenderHints(QPainter::Antialiasing);

	setWindowTitle("Line Segments");
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
				ostrVert << "lines2d.vertices." << vertidx;

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
			for(const Vertex* vertex : m_scene->GetVertexElems())
			{
				QPointF vertexpos = vertex->scenePos();

				std::ostringstream ostrX, ostrY;
				ostrX << "lines2d.vertices." << vertidx << ".<xmlattr>.x";
				ostrY << "lines2d.vertices." << vertidx << ".<xmlattr>.y";

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


	QAction *actionExportGraph = new QAction{"Export Voronoi Graph...", this};
	connect(actionExportGraph, &QAction::triggered, [this]()
	{
		if(QString file = QFileDialog::getSaveFileName(this, "Export DOT", "",
			"DOT Files (*.dot);;All Files (* *.*)"); file!="")
		{
			const auto& graph = m_scene->GetVoroGraph();

			std::ofstream ofstr(file.toStdString());
			print_graph(graph, ofstr);
			ofstr << std::endl;
		}
	});


	QAction *actionQuit = new QAction{"Exit", this};
	connect(actionQuit, &QAction::triggered, [this]() { this->close(); });


	QAction *actionVoronoiRegions = new QAction{"Voronoi Regions", this};
	actionVoronoiRegions->setCheckable(true);
	actionVoronoiRegions->setChecked(true);
	connect(actionVoronoiRegions, &QAction::toggled, [this](bool b)
	{
		m_scene->SetCalcVoro(b);
		m_scene->UpdateVoro();
	});

	QAction *actionVoronoiHelpers = new QAction{"Voronoi Helper Lines", this};
	actionVoronoiHelpers->setCheckable(true);
	actionVoronoiHelpers->setChecked(false);
	connect(actionVoronoiHelpers, &QAction::toggled, [this](bool b)
	{
		m_scene->SetCalcVoroHelpers(b);
		m_scene->UpdateVoro();
	});

	QAction *actionVoronoiRegionsPixel = new QAction{"Voronoi Regions (Pixel-Wise)", this};
	connect(actionVoronoiRegionsPixel, &QAction::triggered, [this]()
	{
		m_scene->UpdateVoroImage(m_view->viewportTransform());
	});


	// intersection backend
	QAction *actionIntersDirect = new QAction{"Direct", this};
	actionIntersDirect->setCheckable(true);
	actionIntersDirect->setChecked(false);
	connect(actionIntersDirect, &QAction::toggled, [this]()
	{ m_scene->SetIntersectionCalculationMethod(IntersectionCalculationMethod::DIRECT); });

	QAction *actionIntersSweep = new QAction{"Sweep", this};
	actionIntersSweep->setCheckable(true);
	actionIntersSweep->setChecked(true);
	connect(actionIntersSweep, &QAction::toggled, [this]()
	{ m_scene->SetIntersectionCalculationMethod(IntersectionCalculationMethod::SWEEP); });


	// voronoi backend
	QAction *actionVoroBoost = new QAction{"Boost", this};
	actionVoroBoost->setCheckable(true);
	actionVoroBoost->setChecked(false);
	connect(actionVoroBoost, &QAction::toggled, [this]()
	{ m_scene->SetVoronoiCalculationMethod(VoronoiCalculationMethod::BOOST); });

	QAction *actionVoroCGAL = new QAction{"CGAL", this};
	actionVoroCGAL->setCheckable(true);
	actionVoroCGAL->setChecked(true);
	connect(actionVoroCGAL, &QAction::toggled, [this]()
	{ m_scene->SetVoronoiCalculationMethod(VoronoiCalculationMethod::CGAL); });


	// groups
	QActionGroup *groupInters = new QActionGroup{this};
	groupInters->addAction(actionIntersDirect);
	groupInters->addAction(actionIntersSweep);

	QActionGroup *groupVoro = new QActionGroup{this};
	groupVoro->addAction(actionVoroBoost);
	groupVoro->addAction(actionVoroCGAL);


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
	menuFile->addAction(actionExportGraph);
	menuFile->addSeparator();
	menuFile->addAction(actionQuit);

	menuCalc->addAction(actionVoronoiRegions);
	menuCalc->addAction(actionVoronoiHelpers);
	menuCalc->addSeparator();
	menuCalc->addAction(actionVoronoiRegionsPixel);

	menuBack->addAction(actionIntersDirect);
	menuBack->addAction(actionIntersSweep);
	menuBack->addSeparator();
	menuBack->addAction(actionVoroBoost);
	menuBack->addAction(actionVoroCGAL);


	// menu bar
	QMenuBar *menuBar = new QMenuBar{this};
	menuBar->setNativeMenuBar(false);
	menuBar->addMenu(menuFile);
	menuBar->addMenu(menuCalc);
	menuBar->addMenu(menuBack);
	setMenuBar(menuBar);


	// connections
	connect(m_view.get(), &LinesView::SignalMouseCoordinates,
	[this](double x, double y, double vpx, double vpy) -> void
	{
		SetStatusMessage(QString("Scene: x=%1, y=%2, Viewport: x=%3, y=%4.")
			.arg(x, 5).arg(y, 5).arg(vpx, 5).arg(vpy, 5));
	});


	SetStatusMessage("Ready.");
}


void LinesWnd::SetStatusMessage(const QString& msg)
{
	m_statusLabel->setText(msg);
}


void LinesWnd::closeEvent(QCloseEvent *e)
{
	// ------------------------------------------------------------------------
	// save settings
	QSettings settings{this};

	QByteArray geo{this->saveGeometry()}, state{this->saveState()};
	settings.setValue("wnd_geo", geo);
	settings.setValue("wnd_state", state);
	// ------------------------------------------------------------------------

	QMainWindow::closeEvent(e);
}


LinesWnd::~LinesWnd()
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
		app->setApplicationName("lines2d");
		set_locales();

		auto vis = std::make_unique<LinesWnd>();
		vis->show();

		return app->exec();
	}
	catch(const std::exception& ex)
	{
		std::cerr << ex.what() << std::endl;
	}

	return -1;
}
// ----------------------------------------------------------------------------

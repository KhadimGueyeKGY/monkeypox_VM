(window.webpackJsonp=window.webpackJsonp||[]).push([[16],{1179:function(e,t,r){"use strict";r.r(t);var n=r(1),o=r.n(n),a=r(796),i=r(831),l=r(818),u=r(834),c=r(940),s=r(848),f=r(849),d=r(843),b=r(884),m=r(822),h=r(954);function y(e,t){return function(e){if(Array.isArray(e))return e}(e)||function(e,t){var r=null==e?null:"undefined"!=typeof Symbol&&e[Symbol.iterator]||e["@@iterator"];if(null!=r){var n,o,a,i,l=[],u=!0,c=!1;try{if(a=(r=r.call(e)).next,0===t){if(Object(r)!==r)return;u=!1}else for(;!(u=(n=a.call(r)).done)&&(l.push(n.value),l.length!==t);u=!0);}catch(e){c=!0,o=e}finally{try{if(!u&&null!=r.return&&(i=r.return(),Object(i)!==i))return}finally{if(c)throw o}}return l}}(e,t)||p(e,t)||function(){throw new TypeError("Invalid attempt to destructure non-iterable instance.\nIn order to be iterable, non-array objects must have a [Symbol.iterator]() method.")}()}function v(e){return function(e){if(Array.isArray(e))return j(e)}(e)||function(e){if("undefined"!=typeof Symbol&&null!=e[Symbol.iterator]||null!=e["@@iterator"])return Array.from(e)}(e)||p(e)||function(){throw new TypeError("Invalid attempt to spread non-iterable instance.\nIn order to be iterable, non-array objects must have a [Symbol.iterator]() method.")}()}function p(e,t){if(e){if("string"==typeof e)return j(e,t);var r=Object.prototype.toString.call(e).slice(8,-1);return"Object"===r&&e.constructor&&(r=e.constructor.name),"Map"===r||"Set"===r?Array.from(e):"Arguments"===r||/^(?:Ui|I)nt(?:8|16|32)(?:Clamped)?Array$/.test(r)?j(e,t):void 0}}function j(e,t){(null==t||t>e.length)&&(t=e.length);for(var r=0,n=new Array(t);r<t;r++)n[r]=e[r];return n}var O=function(e){var t=e.hoverData;if(null===t)return null;var r=t.hoverTitle,n=t.mouseX,a=t.mouseY,i=t.containerId,l=t.data,u={position:"absolute",minWidth:200,padding:"5px",borderRadius:10,zIndex:1e3,pointerEvents:"none",fontFamily:m.e.panel.fontFamily,fontSize:14,lineHeight:1,fontWeight:m.e.panel.fontWeight,color:m.e.panel.color,backgroundColor:m.e.panel.backgroundColor,wordWrap:"break-word",wordBreak:"break-word"},c=document.getElementById(i),s=c.getBoundingClientRect();u.maxWidth=.5*s.width;var f=a-s.top+c.scrollTop,d=n-s.left;return d<.5*s.width?u.left=d+5:u.right=s.width-d+5,f-c.scrollTop<.5*s.height?u.top=f+5:u.bottom=s.height-f+5,o.a.createElement("div",{style:u},o.a.createElement("div",{className:"tooltip",style:m.e.tooltip},o.a.createElement("div",{style:m.e.tooltipHeading},r),v(l.entries()).map((function(e){var t=y(e,2),r=t[0],n=t[1];return o.a.createElement(h.a,{key:r,name:"".concat(r,":"),value:n})}))))},g=r(941);function w(e){return function(e){if(Array.isArray(e))return A(e)}(e)||function(e){if("undefined"!=typeof Symbol&&null!=e[Symbol.iterator]||null!=e["@@iterator"])return Array.from(e)}(e)||S(e)||function(){throw new TypeError("Invalid attempt to spread non-iterable instance.\nIn order to be iterable, non-array objects must have a [Symbol.iterator]() method.")}()}function E(e,t){return function(e){if(Array.isArray(e))return e}(e)||function(e,t){var r=null==e?null:"undefined"!=typeof Symbol&&e[Symbol.iterator]||e["@@iterator"];if(null!=r){var n,o,a,i,l=[],u=!0,c=!1;try{if(a=(r=r.call(e)).next,0===t){if(Object(r)!==r)return;u=!1}else for(;!(u=(n=a.call(r)).done)&&(l.push(n.value),l.length!==t);u=!0);}catch(e){c=!0,o=e}finally{try{if(!u&&null!=r.return&&(i=r.return(),Object(i)!==i))return}finally{if(c)throw o}}return l}}(e,t)||S(e,t)||function(){throw new TypeError("Invalid attempt to destructure non-iterable instance.\nIn order to be iterable, non-array objects must have a [Symbol.iterator]() method.")}()}function S(e,t){if(e){if("string"==typeof e)return A(e,t);var r=Object.prototype.toString.call(e).slice(8,-1);return"Object"===r&&e.constructor&&(r=e.constructor.name),"Map"===r||"Set"===r?Array.from(e):"Arguments"===r||/^(?:Ui|I)nt(?:8|16|32)(?:Clamped)?Array$/.test(r)?A(e,t):void 0}}function A(e,t){(null==t||t>e.length)&&(t=e.length);for(var r=0,n=new Array(t);r<t;r++)n[r]=e[r];return n}var x=function(e,t,r){var n={};return Object.entries(r).forEach((function(e){var t=E(e,2),r=t[0],o=t[1];n[r]=n[r]||[],o.forEach((function(e,t){e.active&&n[r].push(t)}))})),{activeFilters:n,filteredMeasurements:e.filter((function(e){if(t[e.strain]!==l.b)return!1;for(var r=0,o=Object.entries(n);r<o.length;r++){var a=E(o[r],2),i=a[0],u=a[1];if(u.length>0&&!u.includes(e[i]))return!1}return!0}))}},I=function(e){var t=e.height,r=e.width,l=e.showLegend,f=e.setPanelTitle,d=Object(a.d)((function(e){return function(e){var t=e.tree,r=e.controls.colorScale;return t.nodes.reduce((function(e,n,o){var a=e.treeStrainVisibility,i=e.treeStrainColors;if(!n.hasChildren){a[n.name]=t.visibility[o];var l=Object(u.i)(n,r);if(r.continuous){var s=r.visibleLegendValues.find((function(e){return Object(c.b)(e,n,r)}));void 0!==s&&(l=s)}i[n.name]={attribute:l,color:t.nodeColors[o]}}return e}),{treeStrainVisibility:{},treeStrainColors:{}})}(e)}),i.isEqual),m=d.treeStrainVisibility,h=d.treeStrainColors,y=Object(a.d)((function(e){return e.controls.colorScale.legendValues}),i.isEqual),v=Object(a.d)((function(e){return e.metadata.colorings})),p=Object(a.d)((function(e){return e.controls.colorBy})),j=Object(a.d)((function(e){return e.controls.measurementsGroupBy})),S=Object(a.d)((function(e){return e.controls.measurementsFilters})),A=Object(a.d)((function(e){return e.controls.measurementsDisplay})),I=Object(a.d)((function(e){return e.controls.measurementsShowOverallMean})),M=Object(a.d)((function(e){return e.controls.measurementsShowThreshold})),k=Object(a.d)((function(e){return e.measurements.collectionToDisplay}),i.isEqual),C=k.title,T=k.x_axis_label,F=k.threshold,V=k.fields,B=k.measurements,D=k.groupings,W=Object(n.useRef)(null),q=E(Object(n.useState)(null),2),L=q[0],R=q[1],_=x(B,m,S),U=_.activeFilters,Y=_.filteredMeasurements,z=D.get(j).values,G=z;U[j]&&U[j].length&&(G=U[j]);var H,J,N=Object(g.j)(Y,j,G),P=Object(n.useMemo)((function(){return Object(g.f)(r,B)}),[r,B]),X=Object(n.useMemo)((function(){return Object(g.g)()}),[]),$=(H={xScale:P,yScale:X,x_axis_label:T,threshold:F,groupingOrderedValues:z,groupedMeasurements:N},J=Object(n.useRef)(),Object(i.isEqual)(H,J.current)||(J.current=H),J.current),K=Object(n.useMemo)((function(){return function(e,t,r,n){var o=arguments.length>4&&void 0!==arguments[4]?arguments[4]:null,a=null;if(null!==e){var l=null!==o?"Color by ".concat(Object(u.f)(v,p)," : ").concat(o):null,c=new Map;if("measurement"===t){var s=Object.keys(e).filter((function(e){return V.has(e)})),f=w(V.keys()),d=Object(i.orderBy)(s,(function(e){return f.indexOf(e)}));d.forEach((function(t){c.set(V.get(t).title,e[t])}))}else"mean"===t?(c.set("mean",e.mean.toFixed(2)),c.set("standard deviation",e.standardDeviation?e.standardDeviation.toFixed(2):"N/A")):(console.error('"Unknown data type for hover panel: '.concat(t)),Object.entries(e).forEach((function(e){var t=E(e,2),r=t[0],n=t[1];return c.set(r,n)})));a={hoverTitle:l,mouseX:r,mouseY:n,containerId:g.l,data:c}}R(a)}}),[V,v,p]);Object(n.useEffect)((function(){f("".concat(C||"Measurements"," (grouped by ").concat(V.get(j).title,")"))}),[f,C,V,j]),Object(n.useEffect)((function(){Object(g.d)(W.current),Object(g.i)(W.current,$)}),[$]),Object(n.useEffect)((function(){Object(g.a)(W.current,h),Object(g.e)(W.current,h),Object(g.h)(W.current,$,h,y),Object(g.b)(W.current,K,h)}),[$,h,y,K]),Object(n.useEffect)((function(){Object(g.c)(W.current,A)}),[$,h,y,K,A]),Object(n.useEffect)((function(){Object(g.m)(W.current,"overallMean",I)}),[$,I]),Object(n.useEffect)((function(){Object(g.m)(W.current,"threshold",M)}),[$,M]);return o.a.createElement(o.a.Fragment,null,l&&o.a.createElement(s.a,null,o.a.createElement(b.a,{right:!0,width:r})),o.a.createElement("div",{id:g.l,style:{overflowY:"auto",position:"relative",height:t,width:r}},L&&o.a.createElement(O,{hoverData:L}),o.a.createElement("svg",{id:"d3MeasurementsSVG",width:"100%",ref:W})))};t.default=function(e){var t=e.height,r=e.width,i=e.showLegend,l=Object(a.d)((function(e){return e.measurements.loaded})),u=Object(a.d)((function(e){return e.measurements.error})),c=E(Object(n.useState)("Measurements"),2),s=c[0],b=c[1];return o.a.createElement(d.a,{title:s,titleStyles:{width:r,display:"block",overflow:"hidden",whiteSpace:"nowrap",textOverflow:"ellipsis"}},l&&(u?o.a.createElement(f.a,{style:{height:t,width:r},direction:"column",justifyContent:"center"},o.a.createElement("p",{style:{textAlign:"center"}},u)):o.a.createElement(I,{height:t,width:r,showLegend:i,setPanelTitle:b})))}}}]);